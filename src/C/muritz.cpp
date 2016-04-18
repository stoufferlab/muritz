
// autotools
//#include <config.h>

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

// for option parsing
#include <unistd.h>

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
double nullcost;


void help(){
	cerr << "Incorrect usage. Please RTFM.\n";
	exit(1);
}

int main(int argc, char *argv[])
{
    // relevant parameters for simulated annealing
    void (*printfunc)(void*) = NULL;
    bool pairs = false, randomstart = false;
    double iters_fixed_T = 1.0;
    double t_initial = -1.0;
    double mu_t = 1.001;
    double t_min = 1E-7;
    long degree = 0;
    long cost_function = 2;
    int overlap=2;


    // set the above parameters with command line options
    int opt;
    while((opt = getopt(argc, argv, "vprn:t:c:m:k:l:o:u:")) != -1) {
    	switch (opt) {
    		case 'v':
    			printfunc = &alignment_print;
    			break;
    		case 'p':
    			pairs = true;
    			break;
            case 'r':
                randomstart = true;
                break;
    		case 'n':
    			if(optarg)
    				iters_fixed_T = strtod(optarg, NULL);
    			else
    				help();
    			break;
    		case 't':
    			if(optarg)
    				t_initial = strtod(optarg, NULL);
    			else
    				help();
    			break;
    		case 'c':
    			if(optarg)
    				mu_t = strtod(optarg, NULL);
    			else
    				help();
    			break;
    		case 'm':
    			if(optarg)
    				t_min = strtod(optarg, NULL);
    			else
    				help();
    			break;
    		case 'k':
    			if(optarg)
    				degree = strtoul(optarg, NULL, 0);
    			else
    				help();
    			break;
    		case 'l':
    			if(optarg)
    				cost_function = strtoul(optarg, NULL, 0);
    			else
    				help();
    			break;
    		case 'o':
    			if(optarg)
    				overlap = strtoul(optarg, NULL, 0);
    			else
    				help();
    			break;
    		case 'u':
    			if(optarg)
    				nullcost = strtod(optarg, NULL);
    			else
    				help();
    			break;
    		default: // '?' //
    			help();
    	}
    }

	// set up the random number generator
	gsl_rng_env_setup();
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);

	// read in two files of networks
	read_alignment_data(' ',n1,n2);

  	// set up the alignment between networks
	Alignment * alignment = setup_alignment();
    if(randomstart)
        randomize_alignment(r,alignment);

    // decide on what the node-to-node distance function is

    if(cost_function==0){
        alignment->dfunc = &role_euclidean_distance;
    }else{
        if(cost_function==1){
            alignment->dfunc = &role_correlation;
        }else{
            alignment->dfunc = &role_chisquared;
        }
    }


    // assign simulated annleaing parameters to pass to the function below
    alignment->iters_fixed_T = iters_fixed_T;
    alignment->t_initial = t_initial;
    alignment->mu_t = mu_t;
    alignment->t_min = t_min;
    alignment->degree = degree;

	// set up the simulated annealing parameters
	gsl_siman_params_t params = alignment_params(r,alignment);

	// use simulated annealing to find an optimal alignment
	// print out all of the incremental steps in the the optimization
	gsl_siman_solve(r,
					alignment,
					alignment_energy,
					alignment_step,
					alignment_distance,
					printfunc,
					_copy,
					_copy_construct,
					_destroy,
					sizeof(alignment),
					params);
  	
	// print out the "optimal" alignment
    
    //alignment_print_json(alignment, true, pairs);

    if(overlap!=0 && overlap!=-1 && overlap!=1){

	if(!pairs){
    		cout << "optimal ="; alignment_print(alignment); cout << endl;
	}else{
		cout << "optimal ="; alignment_print_pairs(alignment); cout << endl;
	}
	cout << "energy = " << alignment_energy(alignment) << endl;

    }else{
	overlap_pairs(alignment, pairs, overlap);
	cout << "energy = " << alignment_energy(alignment) << endl;
    }

	// free allocated memory
	alignment_free(alignment);
	gsl_rng_free(r);

	return 0;
}
