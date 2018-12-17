import os
import random
import sys
import tempfile
import subprocess
import operator

from pymfinder import motif_roles
from optionparser import parse_cl_options, UNIPARTITE_MOTIF_SIZE, BIPARTITE_MOTIF_SIZE
import muritzex


def read_network(filename):
	network = []
	predators=set([])
	prey=set([])
	spe=set([])
	unipartite=False
	inFile = open(filename,'r')
	for line in inFile:
		interaction=line.strip().split()
                spe.add(interaction[0])
                spe.add(interaction[1])
		if unipartite==False and interaction[0] not in prey and interaction[1] not in predators:
			predators.add(interaction[0])
			prey.add(interaction[1])
		else:
			unipartite=True
		network.append(interaction)
	inFile.close()
	return network, unipartite, spe

def muritz_input(net1, net2, net1_roles, net2_roles, pairs):
	all_roles = net1_roles[net1_roles.keys()[0]].keys()

	input = []
	mat = [0 for i in range(5)]
	mat[0] = ('\n'.join([' '.join(link) for link in net1]))
	mat[1] = ('\n'.join([' '.join(map(str,[n]+[net1_roles[n][r] for r in all_roles])) for n in net1_roles]))
	mat[2] = ('\n'.join([' '.join(link) for link in net2]))
	mat[3] = ('\n'.join([' '.join(map(str,[n]+[net2_roles[n][r] for r in all_roles])) for n in net2_roles]))
	mat[4] = ('\n'.join([' '.join(link) for link in pairs]))
	return(mat)

def class_to_dict(roles): 
    #sort the role vectors
    for n in roles.nodes:
        sorted_list = sorted(roles.nodes[n].roles.items(), key = operator.itemgetter(0))
        roles.nodes[n].roles = {item[0]: item[1] for item in sorted_list}
       
        if roles.weighted:
            #convert from class structure to dict
            class_to_dict = {roles.nodes[n].id:roles.nodes[n].weighted_roles for n in roles.nodes}
        else:
            class_to_dict = {roles.nodes[n].id:roles.nodes[n].roles for n in roles.nodes}

    return class_to_dict

def read_roles(filename, spe):
    inFile=open(filename, "r+")
    lines = [x.strip().split() for x in inFile.readlines()]
    nspe=len(lines)
    inFile.close()
    roles = {x[0]:{idy:float(y) for idy, y in enumerate(x[1:])} for x in lines}

    if nspe!=len(roles.keys()):
        sys.stderr.write("There is something odd about file "+filename+", you should check the name of the nodes and look for repeated names.\n")
        sys.exit()

    if sorted(roles.keys())!=sorted(spe):
        sys.stderr.write("The nodes in the network file do not match the nodes in "+filename+".\n")
        sys.exit()

    if any([len(roles[x])!=len(roles[roles.keys()[0]]) for x in roles.keys()]):
        sys.stderr.write("There is something odd about file "+filename+", all roles should have the same size.\n")
        sys.exit()

    return roles

def muritz(options, args):
    pairs = ""
    pairtype = None

    if(options.fixed_file):
        filename = options.fixed_file 
        pairs, pairtype, sp_ = read_network(filename)

    if(options.roles1!=None and options.roles2==None) or (options.roles1==None and options.roles2!=None):
        sys.stderr.write("If you define the roles, you need to do it for both networks.\n")
        sys.exit()

    net1, net1type, spe1 = read_network(args[0])
    net2, net2type, spe2 = read_network(args[1])

    #Are the networks unipartite or bipartite?
    if options.bipartite:
        if net1type==False and net2type==False:
            unipartite=False
        elif net1type != net2type:
            print "You are comparing a unipartite network with a bipartite one. Both will be considered as unipartite."
            unipartite=True
        else:
            print "Both networks are unipartite, but the bipartite flag was specified. Both will be treated as unipartite."
            unipartite=True
    else:
        unipartite=True

    if options.weighted:
        weighted=True
    else:
        weighted=False
    
    if options.motifsize is None:
        if unipartite:
            motifsize = UNIPARTITE_MOTIF_SIZE
        else:
           motifsize = BIPARTITE_MOTIF_SIZE
    else:
        motifsize = int(options.motifsize)
        if motifsize < 2:
            sys.stderr.write("Motif size must be at least 2.\n")
            sys.exit()
    
    if options.roles1!=None:
        net1_roles=read_roles(options.roles1, spe1)
    else:
        net1_roles = class_to_dict(motif_roles(args[0],motifsize=2, networktype = "unipartite" if unipartite else "bipartite", weighted=weighted, allroles=True))
        for k in range(3,motifsize+1):
            net1_roles_tmp = class_to_dict(motif_roles(args[0],motifsize=k, networktype = "unipartite" if unipartite else "bipartite", weighted=weighted, allroles=True))
            for i in net1_roles:
                net1_roles[i].update(net1_roles_tmp[i])
    
    if options.roles2!=None:
        net2_roles=read_roles(options.roles2, spe2)
    else:
        net2_roles = class_to_dict(motif_roles(args[1],motifsize=2, networktype = "unipartite" if unipartite else "bipartite", weighted=weighted, allroles=True))
        for k in range(3,motifsize+1):
            net2_roles_tmp = class_to_dict(motif_roles(args[1],motifsize=k, networktype = "unipartite" if unipartite else "bipartite", weighted=weighted, allroles=True))
            for i in net2_roles:
                net2_roles[i].update(net2_roles_tmp[i])
    
    if len(net1_roles[net1_roles.keys()[0]])!=len(net2_roles[net2_roles.keys()[0]]):
        sys.stderr.write("The roles for the nodes of each network need to have the same size.\n")
        sys.exit()
    
    muritz_in = muritz_input(net1, net2, net1_roles, net2_roles, pairs)
    
    # get a random seed
    rnd_seed = random.randint(0,sys.maxint)
    os.environ['GSL_RNG_SEED'] = str(rnd_seed)

    # do we want to start from a random initial alignment?
    if options.randomstart:
            rflag = "-r"
    else:
            rflag = ""

    # do we want to track the output from muritz?
    if options.verbose:
            vflag = "-v"
    else:
            vflag = ""

    # do we want to print pairwise energies?the output from muritz?
    if options.pairs:
            pflag = "-p"
    else:
            pflag = ""

    # Did you say bipartite?
    if unipartite:
            bflag = ""
    else:
            bflag = "-b"
    
    # call the muritz alignment code
    #command = "GSL_RNG_SEED=%s muritz.x -n %s -t %s -c %s -m %s -k %s -l %s -o %s -u %s %s %s %s" % (rnd_seed,
    command = "-n %s -t %s -c %s -m %s -k %s -l %s -o %s -u %s -e %s -a %s %s %s %s %s" % (
            options.iterations,
            options.tinitial,
            options.cooling,
            options.tminimum,
            options.degree,
            options.cost_function,
            options.overlap,
            options.nullcost,
            options.endcounter,
            options.acceptmin,
            rflag,
            pflag,
            vflag,
            bflag)
    
    muritzex.muritz(command, muritz_in[0], muritz_in[1], muritz_in[2], muritz_in[3], muritz_in[4])
#	# print out and store the muritz stdout line by line as it comes

#	for line in iter(process.stdout.readline, ''):
#		print(line),

#	process.wait()
		
#	# get rid of the GSL seed info and any empty lines from stderr
#	muritz_stderr = [i for i in process.stderr.readlines() if "GSL_RNG_SEED" not in i and i != '']
#	if muritz_stderr:
#		sys.stderr.write(''.join(muritz_stderr))



	# print out the optimization output
	#print(muritz_stdout)
	#muritz_out.seek(0)
	#print(muritz_out.readlines())
	#for line in muritz_out:
	#	print(line)
	#sys.stdout.write(muritz_stdout)

########################################
########################################
########################################
########################################

def main():
	options, args = parse_cl_options()
	muritz(options, args)

if __name__ == "__main__":
	main()
