##muritz

**muritz** is directed graph alignment based on node-specific motif-role profiles.

At its core, muritz is a Python wrapper for a C/C++ network-alignment code that uses **[pymfinder](http://github.com/stoufferlab/pymfinder)** for the underlying motif analysis.

## Installation instructions

Installation should be relatively painless, though it might require a bit of elbow grease.

1. You must first install [pymfinder](http://github.com/stoufferlab/pymfinder) since muritz cannot run without it. Anyone with access to muritz should have access to pymfinder...

2. Clone the muritz repository to the location of your greatest desire

		git clone git@github.com:stoufferlab/muritz.git

3. Within the cloned repository, run setup.py to install
	
		python setup.py install
   
    If this doesn't work, try adding '--user' at the end or specifying a different '--prefix=/foo'. Beyond a lack of permissions, the most likely source of consternation is a lack of GNU Standard C++ Library, GSL, and/or GSL BLAS on your machine.

4. Assuming you made it this far, you should in fact be able to run muritz

		muritz test/chains-1.net test/chains-2.net 
		
   If you get an error like 'command not found', you likely installed to somewhere outside of your path. This can be remedied by specifying the path before the executable like
   
   		PATH=$PATH:/foo/bin muritz test/chains1.net test/chains-2.net

5. If you want to see full details of the simulated-annealing based optimization, you can also run muritz with verbose mode turned on

		muritz -v test/chains-1.net test/chains-2.net
		
6. If you'd like to try with larger, more complex networks, I recommend

		muritz src/data/unipartite-1.net src/data/unipartite-1.net
