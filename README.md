##muritz

**muritz** is directed graph alignment based on node-specific motif-role profiles.

At its core, muritz is a Python wrapper for a C/C++ network-alignment code that uses **[pymfinder](http://github.com/stoufferlab/pymfinder)** for the underlying motif analysis.

## Installation instructions

Installation is no where close to straightforward and hence not appropriate for the faint of heart. It requires patience, a long (but shortening) series of commands, and maybe a bit of elbow grease. If you've had enough with the idioms and platitudes, but are still there, please proceed as follows.

1. You must first install [pymfinder](http://github.com/stoufferlab/pymfinder) since muritz cannot run without it.

2. Clone the muritz repository. (If you don't know how to do this already, please check github's [Help](https://help.github.com/).)

3. Within the cloned repository, navigate to the 'muritz/src/C/' directory to compile the C/C++ code by running
	
		./autogen.sh
		./configure
		make

4. Assuming you made it this far, you should in fact be able to head back to root directory ('muritz/') and confirm that everything works. Try muritz out first by aligning some simple network motifs (since that will take less time)

		./muritz src/data/motif-1.net src/data/motif-2.net

5. If you want to see full details of the optimization, you can also run muritz with verbose mode

		./muritz -v src/data/motif-1.net src/data/motif-2.net
		
6. If you'd like to try with larger, more complex networks, I can recommend

		./muritz src/data/unipartite-1.net src/data/unipartite-1.net
