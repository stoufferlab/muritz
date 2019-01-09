#!/usr/bin/python

###########################################

# system modules
from optparse import OptionParser
import sys

###########################################

UNIPARTITE_MOTIF_SIZE = 3
BIPARTITE_MOTIF_SIZE = 6

###########################################

def parse_cl_options():
    usage = "usage: %prog [OPTION]... [-x ROLES1_FILE -y ROLES2_FILE] NETWORK_FILE NETWORK_FILE"
    #usage = "usage: %prog [OPTION] FIRST_NETWORK_FILE SECOND_NETWORK_FILE"
    parser = OptionParser(usage)
    
    parser.add_option("-k", "--degree",
                      action="store", dest="degree", type="int",
                      help="degree of alignment to conduct [default: %default]",
                      default=0,
                     )
    parser.add_option("-l", "--cost_function",
                      action="store", dest="cost_function", type="int",
                      help="Euclidean distance (0), Pearson's correlation coefficient (1) or Chi-squared test (2) [default: %default]",
                      default=1,
                     )
    parser.add_option("-t", "--tinitial",
                      action="store", dest="tinitial", type="float",
                      help="initial temperature for simulated annealing (-1 to automatically select) [default: %default]",
                      default=-1,
                     )
    parser.add_option("-m", "--tminimum",
                      action="store", dest="tminimum", type="float",
                      help="minimum temperature for simulated annealing (0 to terminate only on end counter) [default: %default]",
                      default=1E-7,
                     )
    parser.add_option("-c", "--cooling",
                      action="store", dest="cooling", type="float",
                      help="damping factor for temperature [default: %default]",
                      default=1.001,
                     )
    parser.add_option("-n", "--iterations",
                      action="store", dest="iterations", type="float",
                      help="proportion of possible swaps to perform at each temperature [default: %default]",
                      default=1,
                     )
    parser.add_option("-r", "--randomstart",
                      action="store_true", dest="randomstart",
                      help="start from a random alignment [default: %default]",
                      default=False,
                     )
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="print out incremental alignment data [default: %default]",
                      default=False,
                     )
    parser.add_option("-p", "--pairs",
                      action="store_true", dest="pairs",
                      help="print out pairwise alignment data [default: %default]",
                      default=False,
                     )
    parser.add_option("-o", "--overlap",
                      action="store", dest="overlap", type="int",
                      help="calculate the overlap between networks [default: %default]",
                      default=2,
                     )
    parser.add_option("-u", "--nullcost",
                      action="store", dest="nullcost", type="float",
                      help="contribution of non-aligned nodes to the cost function [default: %default]",
                      default=1,
                     )
    parser.add_option("-b", "--bipartite",
                      action="store_true", dest="bipartite",
                      help="consider the networks as bipartite [default: %default]",
                      default=False,
                     )
    
    parser.add_option("-f", "--fixed",
                      action="store", dest="fixed_file", type="string",
                      help="file with fixed pairs [default: %default]",
                      default=False,
                     )
    
    parser.add_option("-w", "--weighted",
                      action="store_true", dest="weighted",
                      help="use interaction strengths to weight alignment [default: %default]",
                      default=False,
                     )
    
    parser.add_option("-e", "--endcounter",
                      action="store", dest="endcounter",
                      help="number of temperatures with very low acceptance rate since the last new best energy before terminating (0 to terminate only at min temperature) [default: %default]",
                      default=100,
                     )
    parser.add_option("-a", "--acceptmin",
                      action="store", dest="acceptmin",
                      help="minimum proportion of steps at a temperature, below which the end counter is incremented [default: %default]",
                      default=0.0,
                     )
    
    parser.add_option("-s", "--motif-size",
                      action="store", dest="motifsize",
                      help="maximum size of motifs to use for roles [default: {} if unipartite, {} if bipartite]".format(UNIPARTITE_MOTIF_SIZE, BIPARTITE_MOTIF_SIZE),
                      default=None,
                     )
    
    parser.add_option("-x", "--first-roles", dest="roles1_file", help="read role data for first network from ROLES1_FILE",)
    parser.add_option("-y", "--second-roles", dest="roles2_file", help="read role data for second network from ROLES2_FILE",)
    parser.set_defaults(roles1=None, roles2=None)
    
    (options, args) = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        sys.exit()
    else:
        return (options, args)
