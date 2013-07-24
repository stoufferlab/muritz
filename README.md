##muritz

**muritz** is directed graph alignment based on node-specific motif-role profiles.

At its core, muritz is a Python wrapper for a C/C++ network alignment code, but that also uses **[pymfinder](http://github.com/stoufferlab/pymfinder) for the underlying motif analyis.

## Installation instructions

Installation is no where close to straightforward and hence not appropriate for the faint of heart. It requires patience, a long (but shortening) series of commands, and maybe a bit of elbow grease. If you've had enough with the idioms and platitudes, but are still there, please proceed as follows.

1. You must first install [pymfinder](http://github.com/stoufferlab/pymfinder) since muritz cannot run without it. Installation instructions for pymfinder *are* complete and should work without any hiccups.

2. Clone the muritz repository. If you don't know how to do this, please check github's [Help](https://help.github.com/).

3. Within the clone repository, navigate to the 'muritz/src/C/' directory to compile the C/C++ code by running
	
	./autogen.sh
	./configure
	make

4. Assuming you made it this far, you should in fact be able to head back to root directory ('muritz/') and confirm that everything works. Try it out by aligning some network motifs first (since that will take less time)

	./muritz src/data/motif-1.net src/data/motif-2.net

If you want to see full details of the optimization, you can also run muritz with verbose mode
If you receive an error about 'Permission denied' or something similar, you most likely don't have permission to install pymfinder in the global Python site-packages or dist-packages directory. In that case, you can install it locally by adding the `--user` option

	./muritz -v src/data/motif-1.net src/data/motif-2.net
