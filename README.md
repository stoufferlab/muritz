# muritz [![](https://badgen.net/badge/DOI/10.1038%2Fs41467-018-05056-0/red)](https://doi.org/10.1038/s41467-018-05056-0)

**muritz** is directed graph alignment based on node-specific motif-role profiles.

At its core, muritz is a Python wrapper for a C/C++ network-alignment code that uses **[pymfinder](http://github.com/stoufferlab/pymfinder)** for the underlying motif analysis. This software was written and tested in Linux machines, but it has been used in both Windows and mac OS.

## Installation instructions

Installation should be relatively painless, though it might require a bit of elbow grease.

1. You must first install [pymfinder](http://github.com/stoufferlab/pymfinder) since muritz cannot run without it. Anyone with access to muritz should have access to pymfinder...

2. Clone the muritz repository to the location of your greatest desire

		git clone git@github.com:stoufferlab/muritz.git

3. Within the cloned repository, run setup.py to install
	
		python setup.py install
   
   If this doesn't work, try adding '--user' at the end or specifying a different '--prefix=/foo'. Beyond a lack of permissions, the most likely source of consternation is a lack of GNU Standard C++ Library, GSL, and/or GSL BLAS on your machine. The installation process might take a few seconds.

4. Assuming you made it this far, you should in fact be able to run muritz

		muritz test/chains-1.net test/chains-2.net 
		
   If you get an error like 'command not found', you likely installed to somewhere outside of your path. This can be remedied by specifying the path before the executable like
   
   		PATH=$PATH:/foo/bin muritz test/chains-1.net test/chains-2.net


   This should return the following:
    > optimal =  [ (A,A2) (B,B1) (C,C1) (NULL,C2) (NULL,A1) (NULL,B2) ]
    > 
    > Energy = 3
    > 
    > Normalized nodes energy = 3.70074e-17

5. If you want to see full details of the simulated-annealing based optimization, you can also run muritz with verbose mode turned on

		muritz -v test/chains-1.net test/chains-2.net
		
6. If you'd like to try with larger, more complex networks, I recommend

		muritz src/data/unipartite-1.net src/data/unipartite-1.net
    
    The run time is largely dependent on the size of the networks. For small networks, this should only take a few seconds; for larger networks, every alignment might take up to a few minutes.

7. To see all options available, run

		muritz

## How to cite?
To use this software, please make sure you cite both Bramon Mora et. al. (*Identifying a common backbone of interactions underlying food webs from different ecosystems. Nature Communications*, 2018) and Bramon Mora et. al. (*pymfinder: a tool for the motif analysis of binary and quantitative complex networks. bioRxiv*,2018)
