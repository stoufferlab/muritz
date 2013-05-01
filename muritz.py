
import subprocess
import sys

from pymfinder import motif_roles,print_roles

def read_network(filename):
	network = []
	inFile = open(netfile1,'r')
	for line in inFile:
		network.append(line.strip().split())
	inFile.close()
	return network

def muritz_input(netfile1, netfile22):
	net1 = read_network(netfile1)
	net2 = read_network(netfile2)

	net1_roles = motif_roles(netfile1,motifsize=3)
	net2_roles = motif_roles(netfile1,motifsize=3)

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

def muritz(netfile1, netfile2):
	muritz_in = muritz_input(netfile1, netfile2)

	# call the muritz alignment code
	command = "./muritz/src/muritz"
	process = subprocess.Popen(command,
    	                       stdin=subprocess.PIPE,
        	                   stdout=subprocess.PIPE,
							   stderr=subprocess.PIPE,
                	           shell=True)

	muritzout, muritzerr = process.communicate(muritz_in)
	print(muritzout)

########################################
########################################
########################################
########################################

netfile1 = "./data/test.net"
netfile2 = "./data/test2.net"

muritz(netfile1, netfile2)



#process.stdin.write("\n".join([' '.join(map(str,[n]+[net2_roles[n][r] for r in all_roles])) for n in net2_roles]))
#process.stdin.close()#process.stdin.write("/#/##/##/#/\n")

#sys.stdout.write('\n'.join(input))

#process.wait()


#for l in process.stdout:
# 	print l

# print line

