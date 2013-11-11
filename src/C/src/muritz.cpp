
// autotools
#include <config.h>

// c++ header files
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

// gsl header files
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_ieee_utils.h>

// my header files
#include <common.hpp>
#include <alignment.hpp>
#include <network.hpp>
#include <roles.hpp>
#include <simulated_annealing.hpp>

// namespaces
using namespace std;

// the networks are stored as global variables
Network n1;
Network n2;

int main(int argc, char *argv[])
{
	unsigned int i,j,k;

	// set up the random number generator
	gsl_rng_env_setup();
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);

	// read in two files of networks
	read_alignment_data(' ',n1,n2);

  	// set up the alignment between networks
	Alignment * alignment = setup_alignment();

    // decide on what the node-to-node distance function is
    alignment->dfunc = &role_correlation;

	// set up the simulated annealing parameters
	gsl_siman_params_t params = alignment_params(alignment);

	// use simulated annealing to find an optimal alignment
	// print out all of the incremental steps in the the optimization
	if(argc==2 && strcmp(argv[1],"-v")==0)
		gsl_siman_solve(r,
						alignment,
						alignment_energy,
						alignment_step,
						alignment_distance,
						alignment_print,
						_copy,
						_copy_construct,
						_destroy,
						sizeof(alignment),
						params);
	// don't print out anything about the optimization
	else
		gsl_siman_solve(r,
						alignment,
						alignment_energy,
						alignment_step,
						alignment_distance,
						NULL,
						_copy,
						_copy_construct,
						_destroy,
						sizeof(alignment),
						params);
  	
	// print out the "optimal" alignment
	cout << "optimal ="; alignment_print(alignment); cout << endl;

	// free allocated memory
	alignment_free(alignment);
	gsl_rng_free(r);

	return 0;
}
