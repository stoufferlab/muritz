
// autotools
//#include <config.h>

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
#include <network.hpp>
#include <alignment.hpp>
#include <simulated_annealing.hpp>

// namespaces
using namespace std;

// externs
extern Network n1;
extern Network n2;

// globals
bool distance_matrix_def = false; // can we access store of node-to-node distances to speed up some calculations?
gsl_matrix * distance_matrix; // node-to-node distances
double* nulldist1; // distances when unaligned
double* nulldist2; // distances when unaligned

// calculate the role-to-role euclidean distance
double role_euclidean_distance(Role *r1, Role *r2){
	double distance = 0;
	for(unsigned int i=0;i<r1->f.size();++i){
		if(r2->name != "NULL")
			distance += (r1->f[i].frequency - r2->f[i].frequency ) * (r1->f[i].frequency - r2->f[i].frequency);
		else
			distance += (r1->f[i].frequency) * (r1->f[i].frequency);
	}
  	return distance;
}

// calculate the role-to-role correlation coefficient
double role_correlation(Role *r1, Role *r2){
    double r;
	double *f1 = (double*) calloc(r1->f.size(), sizeof(double));
	double *f2 = (double*) calloc(r1->f.size(), sizeof(double));

    if(r1->name == "NULL" || r2->name == "NULL")
        r = 0; // this corresponds to complete lack of correlation
    else{
    	for(unsigned int i=0;i<r1->f.size();++i){
	       	f1[i] = r1->f[i].frequency;
		  	f2[i] = r2->f[i].frequency;
    	}

        r = gsl_stats_correlation(f1, 1,
  		      					  f2, 1,
  			       	    		  r1->f.size());
    }

    return 1 - r;
}

// calculate the role-to-role chi-squared test
double role_chisquared(Role *r1, Role *r2){
	unsigned int i,j,nz_cols=0,total=0,df;
	double expected,chisq;

	// a vector for the row sums		
	int rowsums[2] = {0};
	// initialize a vector for the column sums 
	int *colsums = (int*) calloc(r1->f.size(), sizeof(int));
	for(i=0;i<r1->f.size();++i)
		colsums[i] = 0;

	// calculate the row, column, and total sums
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

	return gsl_cdf_chisq_P(chisq, df);
}

// calculate a full role-to-role distance matrix to speed up the SA
void prepare_distance_matrix(double (*dfunc) (Role*,Role*)){
	unsigned int i,j;
    Role null;
    null.name = "NULL";

    // node-to-node distances
    distance_matrix = gsl_matrix_calloc(n1.nodes.size(),n2.nodes.size());
	for(i=0;i<n1.nodes.size();++i)
		for(j=0;j<n2.nodes.size();++j){
			gsl_matrix_set(distance_matrix, i, j, dfunc(&(n1.roles[i]),&(n2.roles[j])));
        }

    // network 1 node distances when unaligned
    nulldist1 = (double*) calloc(n1.nodes.size(), sizeof(double));
    for(i=0;i<n1.nodes.size();++i)
        nulldist1[i] = dfunc(&(n1.roles[i]), &null);

    // network 2 node distances when unaligned
    nulldist2 = (double*) calloc(n2.nodes.size(), sizeof(double));
    for(i=0;i<n2.nodes.size();++i)
        nulldist2[i] = dfunc(&(n2.roles[i]), &null);

    // in the future we can use this info
    distance_matrix_def = true;
}

// calculate the nth degree neighbor lists to speed up the neighborhood-based SA
void prepare_neighbor_data(unsigned int degree){
    unsigned int i;
    set<Node *> nbrs;
    set<Node *>::iterator nbrs_it;

    // set up the degree-th neighbor data for nodes in network 1
    for(i=0;i<n1.nodes.size();++i){
        // add all prey
        n1.nodes[i]->neighbors[degree] = neighbors(&n1, n1.nodes[i], degree, 1);
        // add all predators
        nbrs = neighbors(&n1, n1.nodes[i], degree, -1);
        for(nbrs_it=nbrs.begin();nbrs_it!=nbrs.end();++nbrs_it)
            n1.nodes[i]->neighbors[degree].insert(*nbrs_it);
    }

    // set up the degree-th neighbor data for nodes in network 1
    for(i=0;i<n2.nodes.size();++i){
        // add all prey
        n2.nodes[i]->neighbors[degree] = neighbors(&n2, n2.nodes[i], degree, 1);
        // add all predators
        nbrs = neighbors(&n2, n2.nodes[i], degree, -1);
        for(nbrs_it=nbrs.begin();nbrs_it!=nbrs.end();++nbrs_it)
            n2.nodes[i]->neighbors[degree].insert(*nbrs_it);
    }
}

// calculate the distance between the roles of two nodes
double node_distance(int i, int j, double (*dfunc) (Role*,Role*)){
    // if we don't the quick ref matrices/vectors defined
    if(!distance_matrix_def)
        prepare_distance_matrix(dfunc);
    
    // use the information in the distance matrix
    if(i != -1 && j != -1)
        return gsl_matrix_get(distance_matrix, i, j);
    else{
        if(i == -1)
            return nulldist2[j];
        else
            return nulldist1[i];
    }
}

// search around up to 'degree' connections away from the aligned nodes and calculate the collective alignment there
double neighbor_distance(Alignment *a, unsigned int m, unsigned int degree){
    int l;
    double d = 0;
    int i,j;
    set<Node *> nbr_i, nbr_j;
    set<Node *>::iterator nbr_it;

    // save as ints to avoid confusion later
    i = a->matches[m].first;
    j = a->matches[m].second;

    // i is not null
    if(i != -1){
        // prepare the lists of neighbors if this is the first time this has been run
        if(n1.nodes[i]->neighbors.count(degree) == 0)
            prepare_neighbor_data(degree);

        // save locally to avoid complications later
        nbr_i = n1.nodes[i]->neighbors[degree];
    }

    // j is not null
    if(j != -1){
        // prepare the lists of neighbors if this is the first time this has been run
        if(n2.nodes[j]->neighbors.count(degree) == 0)
            prepare_neighbor_data(degree);
        
        // save locally to avoid complications later
        nbr_j = n2.nodes[j]->neighbors[degree];
    }

    // let the energizing begin!
    // i is not null
    if(i != -1){
        // j is not null
        if(j != -1){
            // align neighbors using the node with the greatest total number as the baseline
            if(nbr_i.size()>=nbr_j.size()){
                // compute the local alignment for all of i's neighbors
                for(nbr_it=nbr_i.begin(); nbr_it!=nbr_i.end(); ++nbr_it){
                    // who is i's neighbor aligned to?
                    l = a->match1[(*nbr_it)->idx];

                    // if l is not null and is also one of j's neighbors
                    if(l != -1 && nbr_j.count(n2.nodes[l]) != 0)
                        d += node_distance((*nbr_it)->idx, l, a->dfunc);
                    // l is null or is not one of j's neighbors
                    else
                        d += node_distance((*nbr_it)->idx, -1, a->dfunc);
                }
            }else{
                // compute the local alignment for all of j's neighbors
                for(nbr_it=nbr_j.begin(); nbr_it!=nbr_j.end(); ++nbr_it){
                    // who is j's neighbor aligned to?
                    l = a->match2[(*nbr_it)->idx];

                    // if l is not null and is also one of i's neighbors
                    if(l != -1 && nbr_i.count(n1.nodes[l]) != 0)
                        d += node_distance(l, (*nbr_it)->idx, a->dfunc);
                    // l is null or is not one of j's neighbors
                    else
                        d += node_distance(-1, (*nbr_it)->idx, a->dfunc);
                }
            }
        }
        // j is null
        else{
            // all neighbors of i are treated as unaligned
            for(nbr_it=nbr_i.begin(); nbr_it!=nbr_i.end(); ++nbr_it){
                d += node_distance((*nbr_it)->idx, -1, a->dfunc);
            }
        }
    }
    // i is null
    else{
        // j is not null
        if(j != -1){
            // have we already calculated j's list of degree-th neighbors?
            if(n2.nodes[j]->neighbors.count(degree) == 0)
                prepare_neighbor_data(degree);

            // save within a local pointer to avoid complications later
            nbr_j = n2.nodes[j]->neighbors[degree];

            // all neighbor nodes are treated as unaligned (i is null)
            for(nbr_it=nbr_j.begin(); nbr_it!=nbr_j.end(); ++nbr_it){
                d += node_distance(-1, (*nbr_it)->idx, a->dfunc);
            }
        }
        // j is null
        else{
            d = 0; // for goodness sake...
        }
    }

    return d;
}

// calculate the weighted distance between two nodes based on the overall alignment
// TODO: set this function up so that we can give it a vector of weights across different "neighborness"
double distance(Alignment *a, unsigned int i){
    double d;
    if(a->degree == 0)
        d = node_distance(a->matches[i].first, a->matches[i].second, a->dfunc);
    else
        d = neighbor_distance(a,i,a->degree);
    return d;
}

// set up the SA parameter values
// TODO: this should be made far more refined by actually using the data to inform the SA
gsl_siman_params_t alignment_params(const gsl_rng * r, void *xp){
    // typecast the alignment object
	Alignment * a = (Alignment *) xp;

    // SA parameter struct
	gsl_siman_params_t params;
	
	// max step size in random walk
	params.step_size = 0.0;
    // number of attempts before stepping
    params.n_tries = 2.0;
    // Boltzmann constant
    params.k = 1.0;
        	
    // number of iterations at each temperature
    params.iters_fixed_T = int((a->iters_fixed_T) * gsl_pow_2(a->matches.size()) + 0.5);
	
    // initial temperature
    if(a->t_initial != -1)
        params.t_initial = a->t_initial;
    else{
        // calculate the average initial change in energy and use it to set the initial temperature
        Alignment * b = setup_alignment();
        _copy(a,b);
        double ae, ae2, de, mean_de, max_de;
        mean_de = 0;
        max_de = 0;
        ae = alignment_energy(b);
        unsigned long shuffles = b->matches.size();
        for(unsigned long i=0;i<shuffles;++i){
            ae2 = ae;
            alignment_step(r,b,0);
            ae = alignment_energy(b);
            de = abs(ae - ae2);
            mean_de += de;
            max_de = max(max_de, de);
        }
        mean_de = de/double(shuffles);
        params.t_initial = max_de/0.7;
        alignment_free(b);
    }

    // damping factor for temperature
	params.mu_t = a->mu_t;

    // minimum temperature
	params.t_min = a->t_min;

	return params;
}

// calculate the energy/cost function of an alignment
double alignment_energy(void *xp){
	double E = 0;

	// cast the void parameter as an alignment data type
	Alignment * a = (Alignment *) xp;

    // sum the cost function across all paired and unpaired nodes
	for(unsigned int i=0;i<a->matches.size();++i){
		E += distance(a, i);
	}

	return E;
}

/* make a move in the alignment space */
void alignment_step(const gsl_rng * r, void *xp, double step_size){
	// prevent warnings about unused parameter
	double dummy = step_size;
    if(dummy != step_size){
        cerr << "Who framed Roger Rabbit?\n" << endl;
        exit(1);
    }

	// case the alignment as an alignment
	Alignment * a = (Alignment *) xp;

	// pick the pairs to swap
	unsigned int p1 = gsl_rng_uniform_int(r,a->matches.size());
	unsigned int p2 = gsl_rng_uniform_int(r,a->matches.size());

	// swap the indices for net2 within the core alignment object
	unsigned int tmp = a->matches[p1].second;
	a->matches[p1].second = a->matches[p2].second;
	a->matches[p2].second = tmp;

    // swap the indices in the first cheater alignment object
    if(a->matches[p1].first != -1)
        a->match1[a->matches[p1].first] = a->matches[p1].second;
    if(a->matches[p2].first != -1)
        a->match1[a->matches[p2].first] = a->matches[p2].second;

    // swap the indices in the second cheater alignment object
    if(a->matches[p1].second != -1)
        a->match2[a->matches[p1].second] = a->matches[p1].first;
    if(a->matches[p2].second != -1)
        a->match2[a->matches[p2].second] = a->matches[p2].first;
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
	unsigned int i;
	int j, k;
	//Role r1, r2;

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

// print out an alignment with pairwise energies
void alignment_print_pairs(void *xp){
	Alignment * a = (Alignment *) xp;
	unsigned int i;
	int j, k;
	//Role r1, r2;

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
            cout << ":";
            cout << distance(a, i);
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

		if(a1->matches[i].first != -1){
			a2->match1[a1->matches[i].first] = a1->matches[i].second;
        }
		if(a1->matches[i].second != -1){
			a2->match2[a1->matches[i].second] = a1->matches[i].first;
        }
    }

    a2->dfunc = a1->dfunc;
    a2->iters_fixed_T = a1->iters_fixed_T;
    a2->t_initial = a1->t_initial;
    a2->mu_t = a1->mu_t;
    a2->t_min = a1->t_min;
    a2->degree = a1->degree;
}

// copy constructor for an alignment
void * _copy_construct(void *xp){
	Alignment *a1 = (Alignment *) xp;
	Alignment *a2 = alignment_alloc(a1->match1.size(),a1->match2.size());
	_copy(a1, a2);
	return a2;
}

// destructor for an alignment
void _destroy(void *xp){
	alignment_free((Alignment *) xp);
}

