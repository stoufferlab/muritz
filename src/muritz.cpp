
// autotools
#include <config.h>

// c++ header files
#include <cstdlib>
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
#include <read_roles.hpp>
#include <simulated_annealing.hpp>

// namespaces
using namespace std;

#define N_TRIES 200             /* how many points do we try before stepping */
#define ITERS_FIXED_T 20      /* how many iterations for each T? */
#define STEP_SIZE 1.0           /* max step size in random walk */
#define K 1.0                   /* Boltzmann constant */
#define T_INITIAL 5000.0        /* initial temperature */
#define MU_T 1.002              /* damping factor for temperature */
#define T_MIN 5.0e-1

gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
                             K, T_INITIAL, MU_T, T_MIN};

int main(int argc, char *argv[])
{
	unsigned int i,j,k;

	gsl_rng_env_setup();
	cout << "seed: " << gsl_rng_default_seed << endl;
	gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);

  	Alignment a_initial;

  	Network n1 = read_roles(argv[1],' ');
  	Network n2 = read_roles(argv[2],' ');

  	if(n1.roles.size() > n2.roles.size()){
  		a_initial.net1 = &n1;
  		a_initial.net2 = &n2;	
  	}else{
  		a_initial.net1 = &n2;
  		a_initial.net2 = &n1;
  	}
  	
  	for(i=0;i<a_initial.net1->roles.size();++i){
  		a_initial.matches.push_back(pair<int,int>(i,-1));
  		a_initial.match1.push_back(i);
  	}

  	for(i=0;i<a_initial.net2->roles.size();++i){
  		a_initial.matches.push_back(pair<int,int>(-1,i));
  	}

	gsl_siman_solve(r, &a_initial,
					alignment_energy,
					alignment_step,
					alignment_distance,
					alignment_print,
					NULL,
					NULL,
					NULL,
					sizeof(a_initial),
					params);
  	
	cout << "best alignment:\n";
	alignment_print(&a_initial);
	cout << alignment_energy(&a_initial) << endl;

	gsl_rng_free (r);
	return 0;
}
