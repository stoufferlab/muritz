
// autotools
#include <config.h>

// c++ header files
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

// gsl header files
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

// local includes
#include <common.hpp>
#include <alignment.hpp>
#include <simulated_annealing.hpp>

// namespaces
using namespace std;

// externs
extern Network n1;
extern Network n2;

// keep some parameters as globals (whether I like it or not)
gsl_matrix * distance_matrix;

// calculate the distance between two roles
double role_euclidean_distance(Role *r1, Role *r2){
	double distance = 0;
	for(int i=0;i<r1->f.size();++i){
		if(r2->name != "NULL")
			distance += (r1->f[i].frequency - r2->f[i].frequency ) * (r1->f[i].frequency - r2->f[i].frequency);
		else
			distance += (r1->f[i].frequency) * (r1->f[i].frequency);
	}
  	return distance;
}

// calculate the correlation between two roles
double role_correlation(Role *r1, Role *r2){
	double f1[r1->f.size()];
	double f2[r1->f.size()];

	for(int i=0;i<r1->f.size();++i){
		f1[i] = r1->f[i].frequency;

		if(r2->name != "NULL")
			f2[i] = r2->f[i].frequency;
		else
			f2[i] = 0;
	}

  	return gsl_stats_correlation(f1, 1,
  								 f2, 1,
  								 r1->f.size());
}

// role to role comparison based on chi-squared statistic
void role_chisquared(Role *r1, Role *r2, double& chisq, int& df){
	int i,j,nz_cols=0,total=0;
	double expected;

	// a vector for the row sums		
	int rowsums[2] = {0};
	// initialize a vector for the column sums 
	int colsums[r1->f.size()];
	for(i=0;i<r1->f.size();++i)
		colsums[i] = 0;

	// a counter for total number of "observations"
	total = 0;
	for(i=0;i<r1->f.size();++i){
		j = r1->f[i].frequency;

		total += j;
		rowsums[0] += j;
		colsums[i] += j;

		if(r2->name != "NULL"){
			j = r2->f[i].frequency;

			total += j;
			rowsums[1] += j;
			colsums[i] += j;
		}	
	}

	// sum the chisquared statistic over columns (rows are hardcoded below)
	chisq = 0;
	for(i=0;i<r1->f.size();++i){
		if(colsums[i] != 0){
			// a column that contributes to the total possible degrees of freedom
			++nz_cols;

			// expected and chisquared contribution for 0,i
			expected = rowsums[0] * colsums[i] / float(total);
			chisq += gsl_pow_2(r1->f[i].frequency - expected) / float(expected);

			if(r2->name != "NULL"){
				// expected and chisquared contribution for 1,i
				expected = rowsums[1] * colsums[i] / float(total);
				chisq += gsl_pow_2(r2->f[i].frequency - expected) / float(expected);
			}
		}
	}

	// calculate the degrees of freedom for the chisquared test 
	// the final values depends on the number of non-zero columns
	df = (nz_cols-1) * (2-1);

	return;
}

// calculate the distance between two roles
double role_distance(Role *r1, Role *r2){
	if(r2->name == "NULL"){
		return 1;
	}
	else{
		double chisq; int df;
		role_chisquared(r1, r2, chisq, df);
		return gsl_cdf_chisq_P(chisq, df);
	}
}

// set up the SA parameter values
// TODO: this can be made far more elegant and refined
gsl_siman_params_t alignment_params(void *xp){
	Alignment * a = (Alignment *) xp;
	
	int N_TRIES = 2.0;             							/* how many points do we try before stepping */
	int ITERS_FIXED_T = gsl_pow_2(a->matches.size());       /* how many iterations for each T? */
	double STEP_SIZE = 0.0;        							/* max step size in random walk */
	double K = 1.0;                							/* Boltzmann constant */
	double T_INITIAL = 1/0.7;                               /* initial temperature */
	double MU_T = 1.001;                                    /* damping factor for temperature */
	double T_MIN = 1.0e-7;							        /* minimum temperature */

	gsl_siman_params_t params;
	
	params.n_tries = N_TRIES;
	params.iters_fixed_T = ITERS_FIXED_T;
	params.step_size = STEP_SIZE;
	
	params.k = K;
	params.t_initial = T_INITIAL;
	params.mu_t = MU_T;
	params.t_min = T_MIN;

	return params;
}

// calculate a full role-to-role distance matrix to speed up the SA
void prepare_distance_matrix(void){
	unsigned int i,j;

	distance_matrix = gsl_matrix_calloc(n1.nodes.size(),n2.nodes.size());

	for(i=0;i<n1.nodes.size();++i)
		for(j=0;j<n2.nodes.size();++j)
			gsl_matrix_set(distance_matrix, i, j, role_distance(&(n1.roles[i]), &(n2.roles[j])));
}

// calculate the energy/cost function of an alignment
double alignment_energy(void *xp){
	Alignment * a = (Alignment *) xp;
	unsigned int i, j, k;
	double E;
	Role null;
	null.name = "NULL";

	E = 0;

	for(i=0;i<a->matches.size();++i){
		j = a->matches[i].first;
		k = a->matches[i].second;

		if(j != -1){
			if(k != -1)
				E += gsl_matrix_get(distance_matrix, j, k);
			else
				E += role_distance(&(n1.roles[j]),&null);
		}else
			if(k != -1)
				E += role_distance(&(n2.roles[k]),&null);
	}

	return E;
}

/* make a move in the alignment space */
void alignment_step(const gsl_rng * r, void *xp, double step_size){
	int i,j;

	// prevent warnings about unused parameter
	step_size = 0 ; 

	// case the alignment as an alignment
	Alignment * a = (Alignment *) xp;

	// pick the pairs to swap
	unsigned int p1 = gsl_rng_uniform_int(r,a->matches.size());
	unsigned int p2 = gsl_rng_uniform_int(r,a->matches.size());

	// swap the indices for net2 within the alignment object
	unsigned int tmp = a->matches[p1].second;
	a->matches[p1].second = a->matches[p2].second;
	a->matches[p2].second = tmp;
}

// calculate the distance between two alignments
double alignment_distance(void *xp, void *yp){
	Alignment *a1 = (Alignment *) xp, *a2 = (Alignment *) yp;
	double distance = 0;
  	for(unsigned int i=0; i<a1->matches.size();++i){
  		// check if each pairwise match is the same
		distance += ((a1->matches[i] == a2->matches[i]) ? 0 : 1);
  	} 
	return distance;
}

// print out an alignment
void alignment_print(void *xp){
	Alignment * a = (Alignment *) xp;
	unsigned int i, j, k;
	Role r1, r2;

	cout << " [";
	for(i=0;i<a->matches.size();++i){
		j = a->matches[i].first;
		k = a->matches[i].second;

		// don't print out NULL matches		
		if(j!=-1 || k!=-1){
			cout << " (";
			if (j != -1)
				cout << n1.roles[j].name;
			else
				cout << "NULL";
			cout << ",";
			if (k != -1)
				cout << n2.roles[k].name;
			else
				cout << "NULL";
			cout <<  ")";
		}
	}
	cout << " ] ";
}

// copy from one alignment to another
void _copy(void *source, void *dest){
	Alignment *a1 = (Alignment *) source, *a2 = (Alignment *) dest;
	for(unsigned int i=0;i<a1->matches.size();++i){
		a2->matches[i] = a1->matches[i];
	}
}

// copy constructor for an alignment
void * _copy_construct(void *xp){
	Alignment *a1 = (Alignment *) xp;
	Alignment *a2 = alignment_alloc(a1->matches.size());
	_copy(a1, a2);
	return a2;
}

// destructor for an alignment
void _destroy(void *xp){
	alignment_free((Alignment *) xp);
}