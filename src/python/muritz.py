import os
import random
import sys
import tempfile
import subprocess
import operator

from pymfinder import motif_roles
from optionparser import parse_cl_options
import muritzex


def read_network(filename):
	network = []
	predators=set([])
	prey=set([])
	unipartite=False
	inFile = open(filename,'r')
	for line in inFile:
		interaction=line.strip().split()
		if unipartite==False and interaction[0] not in prey and interaction[1] not in predators:
			predators.add(interaction[0])
			prey.add(interaction[1])
		else:
			unipartite=True
		network.append(interaction)
	inFile.close()
	return network, unipartite

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
#def fdict(obj):
#    if obj.weighted:
#
#    return obj

def muritz(options, args):
    pairs = ""
    pairtype = None

    if(options.fixed_file):
        filename = options.fixed_file 
        pairs, pairtype = read_network(filename)

    net1, net1type = read_network(args[0])
    net2, net2type = read_network(args[1])

    #Are the networks unipartite or bipartite?
    if options.bipartite:
        if net1type==False and net2type==False:
            unipartite=net1type
        else:
            print "You are comparing a unipartite network with a bipartite one. Both will be considered as unipartite"
            unipartite=True
    else:
    	unipartite=True

    if options.roles1:
        sys.stderr.write("Sorry, the -r option isn't implemented yet.\n")
        sys.exit()
        #net1_roles = read_roles(options.roles1)
    else:
        if unipartite:
            if options.weighted: 
                net1_roles = class_to_dict(motif_roles(args[0],motifsize=2,allroles=True, weighted=True))
                net1_roles2 = class_to_dict(motif_roles(args[0],motifsize=3,allroles=True, weighted=True))
            else: 
                net1_roles = class_to_dict(motif_roles(args[0],motifsize=2,allroles=True))
                net1_roles2 = class_to_dict(motif_roles(args[0],motifsize=3,allroles=True))
                
            for i in net1_roles:
                net1_roles[i].update(net1_roles2[i])
        else:
            net1_roles = motif_roles(args[0],motifsize=2, networktype = "bipartite", allroles=True)
            for k in range(3,5):
                net1_roles2 = motif_roles(args[0],motifsize=k, networktype = "bipartite", allroles=True)
                for i in net1_roles:
                    net1_roles[i].update(net1_roles2[i])


    if options.roles2:
        sys.stderr.write("Sorry, the -s option isn't implemented yet.\n")
        sys.exit()
        #net2_roles = read_roles(options.roles2)
    else:
        if unipartite:
            if options.weighted: 
                net2_roles = class_to_dict(motif_roles(args[1],motifsize=2,allroles=True, weighted=True))
                net2_roles2 = class_to_dict(motif_roles(args[1],motifsize=3,allroles=True, weighted=True))
                print("WEIGHTED!")
                print(net2_roles)
            else: 
                net2_roles = class_to_dict(motif_roles(args[1],motifsize=2,allroles=True))
                net2_roles2 = class_to_dict(motif_roles(args[1],motifsize=3,allroles=True))
                
            for i in net2_roles:
                net2_roles[i].update(net2_roles2[i])
        else:
            net2_roles = motif_roles(args[1],motifsize=2, networktype = "bipartite",allroles=True)
            for k in range(3,5):
                net2_roles2 = motif_roles(args[1],motifsize=k, networktype = "bipartite",allroles=True)
                for i in net2_roles:
                    net2_roles[i].update(net2_roles2[i])
                    
                   
    muritz_in = muritz_input(net1, net2, net1_roles, net2_roles, pairs)
    print(muritz_in)
	
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
    
    # call the muritz alignment code
    #command = "GSL_RNG_SEED=%s muritz.x -n %s -t %s -c %s -m %s -k %s -l %s -o %s -u %s %s %s %s" % (rnd_seed,
    command = "-n %s -t %s -c %s -m %s -k %s -l %s -o %s -u %s %s %s %s" % (
            options.iterations,
            options.tinitial,
            options.cooling,
            options.tminimum,
            options.degree,
            options.cost_function,
            options.overlap,
            options.nullcost,
            rflag,
            pflag,
            vflag)
    
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
