
import os
import random
import subprocess
import sys
import tempfile

from pymfinder import motif_roles,print_roles
from option_parser import parse_cl_options

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

def muritz_input(net1, net2, net1_roles, net2_roles):
	all_roles = net1_roles[net1_roles.keys()[0]].keys()

	input = []

	# write out the links for network 1
	input.append('\n'.join([' '.join(link) for link in net1]))
	input.append('###')

	# write out the roles for network 1
	input.append('\n'.join([' '.join(map(str,[n]+[net1_roles[n][r] for r in all_roles])) for n in net1_roles]))
	input.append('///')

	# write out the links for network 2
	input.append('\n'.join([' '.join(link) for link in net2]))
	input.append('###')

	# write out the roles for network 2
	input.append('\n'.join([' '.join(map(str,[n]+[net2_roles[n][r] for r in all_roles])) for n in net2_roles]))

	return '\n'.join(input)


def muritz(options, args):
	net1, net1type = read_network(args[0])
	net2, net2type = read_network(args[1])

	#Are the networks unipartite or bipartite?
	if net1type==net2type:
		unipartite=net1type
	else:
		print "You are comparing a unipartite network with a bipartite one. Both will be considered as unipartite"
		unipartite=True

	if options.roles1:
		sys.stderr.write("Sorry, the -r option isn't implemented yet.\n")
		sys.exit()
		#net1_roles = read_roles(options.roles1)
	else:
		if unipartite:
			net1_roles = motif_roles(args[0],motifsize=2,)
			net1_roles2 = motif_roles(args[0],motifsize=3,)
			for i in net1_roles:
				net1_roles[i].update(net1_roles2[i])
		else:
			net1_roles = motif_roles(args[0],motifsize=2, networktype = "bipartite",)
			for k in range(3,7):
				net1_roles2 = motif_roles(args[0],motifsize=k, networktype = "bipartite",)
				for i in net1_roles:
					net1_roles[i].update(net1_roles2[i])



	if options.roles2:
		sys.stderr.write("Sorry, the -s option isn't implemented yet.\n")
		sys.exit()
		#net2_roles = read_roles(options.roles2)
	else:
		if unipartite:
			net2_roles = motif_roles(args[1],motifsize=2,)
			net2_roles2 = motif_roles(args[1],motifsize=3,)
			for i in net2_roles:
				net2_roles[i].update(net2_roles2[i])
		else:
			net2_roles = motif_roles(args[1],motifsize=2, networktype = "bipartite",)
			for k in range(3,7):
				net2_roles2 = motif_roles(args[1],motifsize=k, networktype = "bipartite",)
				for i in net2_roles:
					net2_roles[i].update(net2_roles2[i])

	muritz_in = muritz_input(net1, net2, net1_roles, net2_roles)

	# get a random seed
	rnd_seed = random.randint(0,sys.maxint)

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
	command = "GSL_RNG_SEED=%s muritz.x -n %s -t %s -c %s -m %s -k %s -b %s -o %s -u %s %s %s %s" % (rnd_seed,
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
	#muritz_out = tempfile.TemporaryFile()
	process = subprocess.Popen(command,
					bufsize=0,
					stdin=subprocess.PIPE,
					stdout=subprocess.PIPE,
					stderr=subprocess.PIPE,
					shell=True,
					)

	# write the network and role data to muritz
	process.stdin.write(muritz_in)
	process.stdin.close()


	# print out and store the muritz stdout line by line as it comes

	for line in iter(process.stdout.readline, ''):
		print(line),

	process.wait()
		
	# get rid of the GSL seed info and any empty lines from stderr
	muritz_stderr = [i for i in process.stderr.readlines() if "GSL_RNG_SEED" not in i and i != '']
	if muritz_stderr:
		sys.stderr.write(''.join(muritz_stderr))



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
