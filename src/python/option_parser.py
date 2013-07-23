#!/usr/bin/python

###########################################

# system modules
from optparse import OptionParser
import sys

###########################################

def parse_cl_options():
	#usage = "usage: %prog [-r ROLE_FILENAME] [-s ROLE_FILENAME] NETWORK_FILE NETWORK_FILE"
	usage = "usage: %prog FIRST_NETWORK_FILE SECOND_NETWORK_FILE"
	parser = OptionParser(usage)

	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="print out incremental alignment data", default=False)
	#parser.add_option("-r", "--first-roles", dest="roles1", help="read role data for first network from ROLE_FILENAME", )
	#parser.add_option("-s", "--second-roles", dest="roles2", help="read role data for second network from ROLE_FILENAME", )

	parser.set_defaults(roles1=None, roles2=None)

	(options, args) = parser.parse_args()

	if len(args) != 2:
		parser.print_help()
		sys.exit()
	else:
		return (options, args)
