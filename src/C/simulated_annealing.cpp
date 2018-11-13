
// autotools
//#include <config.h>

// c++ header files
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <string>
#include <algorithm>
#include <vector>

// gsl header files
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

// local includes
#include "common.hpp"
#include "network.hpp"
#include "alignment.hpp"
#include "anneal.hpp"
#include "simulated_annealing.hpp"

// namespaces
using namespace std;

// externs
extern Network n1;
extern Network n2;
extern double nullcost;

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
    double rowsums[2] = {0};
    double *f1 = (double*) calloc(r1->f.size(), sizeof(double));
    double *f2 = (double*) calloc(r1->f.size(), sizeof(double));

    if(r1->name == "NULL" || r2->name == "NULL") {
        return nullcost; // this corresponds to complete lack of correlation
    }
    else{
    	for(unsigned int i=0;i<r1->f.size();++i){
	       	f1[i] = r1->f[i].frequency;
            f2[i] = r2->f[i].frequency;
            rowsums[0] += f1[i];
            rowsums[1] += f2[i];
    	}
        if (rowsums[0] == 0 || rowsums[1] == 0){
            return nullcost;
        }else{
            r = gsl_stats_correlation(f1, 1, f2, 1, r1->f.size());
        }
    }
    
    return 1 - r;
}

// calculate the role-to-role chi-squared test
double role_chisquared(Role *r1, Role *r2){
	unsigned int i,j,nz_cols=0,total=0,df;
	double expected,chisq;

	// a vector for the row sums		
	double rowsums[2] = {0};
	// initialize a vector for the column sums 
	int *colsums = (int*) calloc(r1->f.size(), sizeof(int));
	for(i=0;i<r1->f.size();++i)
		colsums[i] = 0;

	if(r2->name == "NULL"){
		return nullcost;
	}else{
		// calculate the row, column, and total sums
		for(i=0;i<r1->f.size();++i){
			j = r1->f[i].frequency;

			total += j;
			rowsums[0] += j;
			colsums[i] += j;

			j = r2->f[i].frequency;

			total += j;
			rowsums[1] += j;
			colsums[i] += j;	
		}

		// sum the chisquared statistic over columns (rows are hardcoded below)
		chisq = 0;
		if (rowsums[0] == 0 || rowsums[1] == 0){
			return nullcost;
		}
		for(i=0;i<r1->f.size();++i){
			if(colsums[i] != 0){
				// a column that contributes to the total possible degrees of freedom
				++nz_cols;

				// expected and chisquared contribution for 0,i
				expected = rowsums[0] * colsums[i] / float(total);
				chisq += gsl_pow_2(r1->f[i].frequency - expected) / float(expected);

				// expected and chisquared contribution for 1,i
				expected = rowsums[1] * colsums[i] / float(total);
				chisq += gsl_pow_2(r2->f[i].frequency - expected) / float(expected);
			}
		}

		// calculate the degrees of freedom for the chisquared test 
		// the final values depends on the number of non-zero columns
		df = (nz_cols-1) * (2-1);

		return gsl_cdf_chisq_P(chisq, df);
	}
}

// calculate a full role-to-role distance matrix to speed up the SA
void prepare_distance_matrix(double (*dfunc) (Role*,Role*)){
	unsigned int i,j;
    Role null;
    null.name = "NULL";

    // node-to-node distances
    distance_matrix = gsl_matrix_calloc(n1.nodes.size(),n2.nodes.size());
	for(i=0;i<n1.nodes.size();++i) {
		for(j=0;j<n2.nodes.size();++j){
			gsl_matrix_set(distance_matrix, i, j, dfunc(&(n1.roles[i]),&(n2.roles[j])));
        }
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

    // set up the degree-th neighbor data for nodes in network 2
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
double node_distance(int i, int j, double (*dfunc) (Role*,Role*)) {
    // if we don't have the quick ref matrices/vectors defined
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

// add the relevant match as one match contributing, and adjust energy accordingly
static void adjust_contributing_matches(Alignment *a, int i, int j, int delta) {
    if(a->matchesContributing.count(make_pair(i, j))) {
        // if it already exists in the map
        a->matchesContributing[make_pair(i, j)] += delta;// add one
    } else {
        // it doesn't yet exist
        a->matchesContributing.insert(make_pair(make_pair(i, j), delta));// insert it with one copy
    }
    
    //adjust energy
    a->energy += node_distance(i, j, a->dfunc);
}

// add a contributing match delta to the map, and adjust proposed energy accordingly
static void adjust_proposed_deltas(Alignment *a, int i, int j, int delta) {
    if(a->proposedContributionDeltas.count(make_pair(i, j))) {
        // if it already exists in the map
        a->proposedContributionDeltas[make_pair(i, j)] += delta;// add delta
    } else {
        // it doesn't yet exist
        a->proposedContributionDeltas.insert(make_pair(make_pair(i, j), delta));// insert it
    }
    
    //adjust energy
    a->proposedEnergy += node_distance(i, j, a->dfunc) * delta;
}

// apply deltas to contributing matches
static void apply_proposed_deltas(Alignment *a) {
    for(map<pair<int, int>, int>::iterator it = a->proposedContributionDeltas.begin(); it != a->proposedContributionDeltas.end(); it++) {
        adjust_contributing_matches(a, (*it).first.first, (*it).first.second, (*it).second);
    }
    a->proposedContributionDeltas.clear();
}

// search around up to 'degree' connections away from the aligned nodes and add matches contributing to the total energy
static void add_matches_contributing(Alignment *a, unsigned int m) {
    unsigned int degree = a->degree;
    int l;
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
    
    // start adding matches
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
                        adjust_contributing_matches(a, (*nbr_it)->idx, l, 1);
                    // l is null or is not one of j's neighbors
                    else
                        adjust_contributing_matches(a, (*nbr_it)->idx, -1, 1);
                }
            }else{
                // compute the local alignment for all of j's neighbors
                for(nbr_it=nbr_j.begin(); nbr_it!=nbr_j.end(); ++nbr_it){
                    // who is j's neighbor aligned to?
                    l = a->match2[(*nbr_it)->idx];
                    
                    // if l is not null and is also one of i's neighbors
                    if(l != -1 && nbr_i.count(n1.nodes[l]) != 0)
                        adjust_contributing_matches(a, l, (*nbr_it)->idx, 1);
                    // l is null or is not one of i's neighbors
                    else
                        adjust_contributing_matches(a, -1, (*nbr_it)->idx, 1);
                }
            }
        }
        // j is null
        else{
            // all neighbors of i are treated as unaligned
            for(nbr_it=nbr_i.begin(); nbr_it!=nbr_i.end(); ++nbr_it){
                adjust_contributing_matches(a, (*nbr_it)->idx, -1, 1);
            }
        }
    }
    // i is null
    else{
        // j is not null
        if(j != -1){
            // all neighbor nodes are treated as unaligned (i is null)
            for(nbr_it=nbr_j.begin(); nbr_it!=nbr_j.end(); ++nbr_it){
                adjust_contributing_matches(a, -1, (*nbr_it)->idx, 1);
            }
        }
        // else i and j are both null, nothing to add
    }
}

// set up the list of matches that contribute to the total energy and find the total energy in the process
// TODO: set this function (and matchesContributing) up so that we can give it a vector of weights across different "neighborness"
static void matches_contributing_setup(Alignment *a){
    a->matchesContributing.clear();// wipe it to prepare to refill it
    a->energy = 0.0;// reset the energy
    
    for(unsigned int i = 0; i < a->matches.size(); i++){
        add_matches_contributing(a,i);
    }
}

// calculate the energy/cost function of an alignment and store it in the alignment
// also sets up matchesContributing, if applicable (if degree != 0)
void alignment_energy_setup(Alignment *a)
{
    if(a->degree == 0) {
        a->energy = 0;
        // sum the cost function across all paired and unpaired nodes
        for(unsigned int i = 0; i < a->matches.size(); i++){
            a->energy += node_distance(a->matches[i].first, a->matches[i].second, a->dfunc);
        }
    } else {// degree != 0
        matches_contributing_setup(a);
    }
    
    /*
    for(int i = 0; i < n1.nodes.size(); i++) {
        cout << "nulldist1[" << n1.roles[i].name << "] = " << nulldist1[i] << endl;
    }
    for(int i = 0; i < n2.nodes.size(); i++) {
        cout << "nulldist2[" << n2.roles[i].name << "] = " << nulldist2[i] << endl;
    }
    for(int i = 0; i < n1.nodes.size(); i++) {
        for(int j = 0; j < n2.nodes.size(); j++) {
            cout << "dist[" << n1.roles[i].name << "][" << n2.roles[j].name << "] = " << gsl_matrix_get(distance_matrix, i, j) << endl;
        }
    }
    //*/
}

// return the energy/cost function of an alignment
double alignment_energy(void *xp)
{
    // cast the void parameter as an alignment data type
    Alignment * a = (Alignment *) xp;
    
    return a->energy;
}

// remove match m, adjusting proposedContributionDeltas and proposed energy
static void propose_remove_match(Alignment *a, int m) {
    unsigned int degree = a->degree;
    int l;
    int i,j;
    set<Node *> nbr_i, nbr_j;
    set<Node *>::iterator nbr_it;
    
    // save as ints to avoid confusion later
    i = a->matches[m].first;
    j = a->matches[m].second;
    
    nbr_i = n1.nodes[i]->neighbors[degree];
    nbr_j = n2.nodes[j]->neighbors[degree];
    
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
                        adjust_proposed_deltas(a, (*nbr_it)->idx, l, -1);
                    // l is null or is not one of j's neighbors
                    else
                        adjust_proposed_deltas(a, (*nbr_it)->idx, -1, -1);
                }
            }else{
                // compute the local alignment for all of j's neighbors
                for(nbr_it=nbr_j.begin(); nbr_it!=nbr_j.end(); ++nbr_it){
                    // who is j's neighbor aligned to?
                    l = a->match2[(*nbr_it)->idx];
                    
                    // if l is not null and is also one of i's neighbors
                    if(l != -1 && nbr_i.count(n1.nodes[l]) != 0)
                        adjust_proposed_deltas(a, l, (*nbr_it)->idx, -1);
                    // l is null or is not one of i's neighbors
                    else
                        adjust_proposed_deltas(a, -1, (*nbr_it)->idx, -1);
                }
            }
        }
        // j is null
        else{
            // all neighbors of i are treated as unaligned
            for(nbr_it=nbr_i.begin(); nbr_it!=nbr_i.end(); ++nbr_it){
                adjust_proposed_deltas(a, (*nbr_it)->idx, -1, -1);
            }
        }
    }
    // i is null
    else{
        // j is not null
        if(j != -1){
            // all neighbor nodes are treated as unaligned (i is null)
            for(nbr_it=nbr_j.begin(); nbr_it!=nbr_j.end(); ++nbr_it){
                adjust_proposed_deltas(a, -1, (*nbr_it)->idx, -1);
            }
        }
        // else i and j are both null, nothing to remove
    }
    
    // Unless it's totally null, also wipe all copies of the match itself
    if(i != -1 || j != -1) {
        adjust_proposed_deltas(a, i, j, -a->matchesContributing[make_pair(i, j)]);
    }
}

// add match m, using proposed matches, and adjusting proposedContributionDeltas and proposed energy
static void propose_add_match(Alignment *a, int m) {
    unsigned int degree = a->degree;
    int l;
    int i,j;
    set<Node *> nbr_i, nbr_j;
    set<Node *>::iterator nbr_it;
    
    //How many times we must add *this* match to the network, because adjacent matches have in essence added it.
    int numMatches = 0;
    
    //Note which nodes have been swapped in which networks
    int net1_s1 = a->matches[a->p1].first;
    int net1_s2 = a->matches[a->p2].second;
    int net2_s1 = a->matches[a->p1].first;
    int net2_s2 = a->matches[a->p2].second;
    
    i = a->matches[m].first;
    j = a->matches[m].second;
    //Adjust i and j for proposed swap. Note that swap is always of net2 nodes.
    if(m == a->p1) {
        j = net2_s2;
    } else if(m == a->p2) {
        j = net2_s1;
    }
    
    nbr_i = n1.nodes[i]->neighbors[degree];
    nbr_j = n2.nodes[j]->neighbors[degree];
    
    // i is not null
    if(i != -1){
        // j is not null
        if(j != -1){
            // align neighbors using the node with the greatest total number as the baseline
            if(nbr_i.size()>=nbr_j.size()){
                // compute the local alignment for all of i's neighbors
                for(nbr_it=nbr_i.begin(); nbr_it!=nbr_i.end(); ++nbr_it){
                    // what is the index of the neighbour we are currently concerned with?
                    int nbr = (*nbr_it)->idx;
                    // who is i's neighbor aligned to, adjusted for proposed swap?
                    if     (nbr == net1_s1) l = net2_s2;
                    else if(nbr == net1_s2) l = net2_s1;
                    else                    l = a->match1[nbr];
                    
                    // if l is not null and is also one of j's neighbors
                    if(l != -1 && nbr_j.count(n2.nodes[l]) != 0) {
                        adjust_proposed_deltas(a, nbr, l, 1);
                        numMatches++;
                    }
                    // l is null or is not one of j's neighbors
                    else
                        adjust_proposed_deltas(a, nbr, -1, 1);
                }
            }else{
                // compute the local alignment for all of j's neighbors
                for(nbr_it=nbr_j.begin(); nbr_it!=nbr_j.end(); ++nbr_it){
                    //what is the index of the neighbour we are currently concerned with?
                    int nbr = (*nbr_it)->idx;
                    // who is j's neighbor aligned to, adjusted for proposed swap?
                    if     (nbr == net2_s1) l = net1_s2;
                    else if(nbr == net2_s2) l = net1_s1;
                    else                    l = a->match2[nbr];
                    
                    // if l is not null and is also one of i's neighbors
                    if(l != -1 && nbr_i.count(n1.nodes[l]) != 0) {
                        adjust_proposed_deltas(a, l, nbr, 1);
                        numMatches++;
                    }
                    // l is null or is not one of i's neighbors
                    else
                        adjust_proposed_deltas(a, -1, nbr, 1);
                }
            }
        }
        // j is null
        else{
            for(nbr_it=nbr_i.begin(); nbr_it!=nbr_i.end(); ++nbr_it){
                // all neighbors of i are treated as unaligned
                adjust_proposed_deltas(a, (*nbr_it)->idx, -1, 1);
                
                // also need to find when *this* match would be added
                // what is the index of the neighbour we are currently concerned with?
                int nbr = (*nbr_it)->idx;
                // who is i's neighbor aligned to, adjusted for proposed swap?
                if     (nbr == net1_s1) l = net2_s2;
                else if(nbr == net1_s2) l = net2_s1;
                else                    l = a->match1[nbr];
                
                // does nbr have more neighbours than its match (net1 wins ties)?
                // its match being null also qualifies
                if(l == -1 ||
                   n1.nodes[nbr]->neighbors[degree].size() >= n2.nodes[l]->neighbors[degree].size()) {
                    numMatches ++;
                }
            }
        }
    }
    // i is null
    else{
        // j is not null
        if(j != -1){
            for(nbr_it=nbr_j.begin(); nbr_it!=nbr_j.end(); ++nbr_it){
                // all neighbor nodes are treated as unaligned (i is null)
                adjust_proposed_deltas(a, -1, (*nbr_it)->idx, 1);
                
                // also need to find when *this* match would be added
                // what is the index of the neighbour we are currently concerned with?
                int nbr = (*nbr_it)->idx;
                // who is j's neighbor aligned to, adjusted for proposed swap?
                if     (nbr == net2_s1) l = net1_s2;
                else if(nbr == net2_s2) l = net1_s1;
                else                    l = a->match2[nbr];
                
                // does nbr have more neighbours than its match (net1 wins ties)?
                // its match being null also qualifies
                if(l == -1 ||
                   n1.nodes[l]->neighbors[degree].size() < n2.nodes[nbr]->neighbors[degree].size()) {
                    numMatches ++;
                }
            }
        }
        // else i and j are both null, nothing to add
    }
    
    // add this match itself the number of times other matches would add it
    adjust_proposed_deltas(a, i, j, numMatches);
}

// calculate the proposed energy from the current energy and the proposed swap
// also updates proposed contributing matches if applicable (if degree != 0)
static void update_proposed_energy(Alignment *a) {
    a->proposedEnergy = a->energy;// Set the current energy to build from.
    if(a->p1 == a->p2) return;// The proposed step does nothing. TODO: Fix stuff so this can't happen in the first place.
    
    if(a->degree != 0) {
        propose_remove_match(a, a->p1);
        propose_remove_match(a, a->p2);
        propose_add_match(a, a->p1);
        propose_add_match(a, a->p2);
    } else {
        a->proposedEnergy -= node_distance(a->matches[a->p1].first, a->matches[a->p1].second, a->dfunc);
        a->proposedEnergy -= node_distance(a->matches[a->p2].first, a->matches[a->p2].second, a->dfunc);
        a->proposedEnergy += node_distance(a->matches[a->p1].first, a->matches[a->p2].second, a->dfunc);
        a->proposedEnergy += node_distance(a->matches[a->p2].first, a->matches[a->p1].second, a->dfunc);
    }
}


/* Propose a move in the alignment space and return the energy of the resulting alignment. */
double alignment_propose_step(void *xp, const gsl_rng *r)
{
    // cast the alignment as an alignment
    Alignment * a = (Alignment *) xp;
    
    // find the probability that the switch is of two nodes within network A, as opposed to B
    static float probA = ((float)a->unfixed_pairs_A.size() / (float)(a->unfixed_pairs_A.size()+a->unfixed_pairs_B.size()));
    
    // pick the pairs to swap
    if(gsl_rng_uniform(r) < probA){
        a->p1 = a->unfixed_pairs_A[gsl_rng_uniform_int(r,a->unfixed_pairs_A.size())];
        a->p2 = a->unfixed_pairs_A[gsl_rng_uniform_int(r,a->unfixed_pairs_A.size())];
    } else {
        a->p1 = a->unfixed_pairs_B[gsl_rng_uniform_int(r,a->unfixed_pairs_B.size())];
        a->p2 = a->unfixed_pairs_B[gsl_rng_uniform_int(r,a->unfixed_pairs_B.size())];
    }
    
    
    // calculate the energy of the new alignment
    update_proposed_energy(a);
    return a->proposedEnergy;
}

/* Commit a proposed step. */
void alignment_commit_step(void *xp)
{
    // cast the alignment as an alignment
    Alignment * a = (Alignment *) xp;
    
    // swap the indices for net2 within the core alignment object
    unsigned int tmp = a->matches[a->p1].second;
    a->matches[a->p1].second = a->matches[a->p2].second;
    a->matches[a->p2].second = tmp;
    
    // swap the indices in the first cheater alignment object
    if(a->matches[a->p1].first != -1)
        a->match1[a->matches[a->p1].first] = a->matches[a->p1].second;
    if(a->matches[a->p2].first != -1)
        a->match1[a->matches[a->p2].first] = a->matches[a->p2].second;
    
    // swap the indices in the second cheater alignment object
    if(a->matches[a->p1].second != -1)
        a->match2[a->matches[a->p1].second] = a->matches[a->p1].first;
    if(a->matches[a->p2].second != -1)
        a->match2[a->matches[a->p2].second] = a->matches[a->p2].first;
    
    a->energy = a->proposedEnergy;
}

/* Propose and immediately make a step. */
void alignment_step(void *xp, const gsl_rng *r)
{
    alignment_propose_step(xp, r);
    alignment_commit_step(xp);
}

// set up the SA parameter values
// TODO: this should be made far more refined by actually using the data to inform the SA
anneal_params_t alignment_params(const gsl_rng *rng,
                                 const Alignment *a,
                                 double initialTemperature,//-1 if we should calculate it now.
                                 double coolingFactor,
                                 double minTemperature,
                                 int stepsPerTemperature) {
    // SA parameter struct
    anneal_params_t params;
    
    // number of iterations at each temperature
    // TODO: Change this next line to:
    // int(stepsPerTemperature * (gsl_pow_2(a->unfixed_pairs_A.size())+gsl_pow_2(a->unfixed_pairs_B.size())) + 0.5); 
    params.stepsPerTemperature = int(stepsPerTemperature * gsl_pow_2(a->unfixed_pairs_A.size()+a->unfixed_pairs_B.size()) + 0.5);
    
    // initial temperature
    if(initialTemperature != -1) {
        params.initialTemperature = initialTemperature;
    } else {
        // calculate the average initial change in energy and use it to set the initial temperature
        Alignment * b = setup_alignment(a->set_pairs);
        _copy(a,b);
        double ae, ae2, de, mean_de, max_de;
        mean_de = 0;
        max_de = 0;
        ae = alignment_energy(b);
        unsigned long shuffles = b->unfixed_pairs_A.size() + b->unfixed_pairs_B.size();
        for(unsigned long i = 0; i < shuffles; i++){
            ae2 = ae;
            alignment_step(b,rng);
            ae = alignment_energy(b);
            de = abs(ae - ae2);
            mean_de += de;
            max_de = max(max_de, de);
        }
        mean_de = de/double(shuffles);
        params.initialTemperature = max_de/0.7;//TODO: Do something about this magic number.
        
        // Ensure at a least some annealing happens. Some seeds never find a change in energy and start temperature at zero otherwise.
        // TODO: Test and improve this. Seed 6315191922508953630 on chains-1 and chains-1 displays this behaviour.
        params.initialTemperature = max(params.initialTemperature, minTemperature*2);
        
        alignment_free(b);
    }
    
    // damping factor for temperature
    params.coolingFactor = coolingFactor;
    
    // minimum temperature
    params.minTemperature = minTemperature;
    
    return params;
}

// print the energy/cost and the normalized energy/cost function of an alignment
// TODO: Sort out the whole printing deal again.
void print_energy(void *xp, int cost_function, long degree){
	/*double E = 0, Epair_nei, Epair_nod, Enorm_nei=0, Enorm_nod=0;
	int j, k, nei1, nei2, norm=0;
	set<Node *> nbr_i;
	
	// cast the void parameter as an alignment data type
	Alignment * a = (Alignment *) xp;
	
	// sum the cost function across all paired and unpaired nodes
	for(unsigned int i=0;i<a->matches.size();++i){
		Epair_nei = distance(a, i);
		Epair_nod = node_distance(a->matches[i].first, a->matches[i].second, a->dfunc);
		E += Epair_nei;
		if (cost_function==1){
			j = a->matches[i].first;
			k = a->matches[i].second;
			if (j != -1 && k != -1){
				Enorm_nod+=Epair_nod;
				norm++;
				if (degree==1){
 					nbr_i = n1.nodes[j]->neighbors[degree];
					nei1=nbr_i.size();
					
 					nbr_i = n2.nodes[k]->neighbors[degree];
					nei2=nbr_i.size();
					
					if (nei1<nei2){
						if (nei1!=0){
							Epair_nei=Epair_nei-nullcost*(nei2-nei1);
							Enorm_nei+=Epair_nei/(double)nei1;
						}else{
							Enorm_nei+=1;
						}
					}else{
						if (nei2!=0){
							Epair_nei=Epair_nei-nullcost*(nei1-nei2);
							Enorm_nei+=Epair_nei/(double)nei2;
						}else{
							Enorm_nei+=1;
						}
					}
				}else{
					Enorm_nei+=Epair_nei;
				}
			}
		}
	}
	
	if (cost_function==1 && nullcost==1){
		Enorm_nei=Enorm_nei/(double)norm;
		Enorm_nod=Enorm_nod/(double)norm;
		cout << "Energy = " << E << endl;
		cout << "Normalized nodes energy = " << Enorm_nod << endl;
		if (degree==1){
			cout << "Normalized neighbours energy = " << Enorm_nei << endl;
		}
	}else{
		cout << "Energy = " << E << endl;
	}*/
	
	printf("Energy = %lf\n", alignment_energy(xp));
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
/*void alignment_print_pairs(void *xp){
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
            cout << ",";
            cout << node_distance(a->matches[i].first, a->matches[i].second, a->dfunc);
            cout <<  ")";
		}
	}
	cout << " ] ";
}*/

/*void overlap_pairs(void *xp, bool pairs, int direction){
	Alignment * a = (Alignment *) xp;
	unsigned int i,m;
	int j, k;
	vector<int> totals(3);
	set<int> v1;
	set<int> v2;
	vector<int> v;
	Node * nbr;


	//Role r1, r2;
	cout << "optimal = [";
	for(i=0;i<a->matches.size();++i){
		j = a->matches[i].first;
		k = a->matches[i].second;
		
		// don't print out NULL matches		
		if(j!=-1 || k!=-1){
			cout << " (";

			if (j != -1 && k != -1){

				cout << n1.roles[j].name;
				cout << ",";
				cout << n2.roles[k].name;

				//find neighbors of species j
				v1.clear();
				nbr = n1.nodes[j];
				if(direction == 0 || direction == 1)
					for(m=0;m<nbr->prey.size();++m)
						v1.insert(nbr->prey[m]->idx);
				if(direction == 0 || direction == -1)
					for(m=0;m<nbr->predators.size();++m)
						v1.insert(nbr->predators[m]->idx);


				//find neighbors of species k
				v2.clear();
				nbr = n2.nodes[k];
				if(direction == 0 || direction == 1)
					for(m=0;m<nbr->prey.size();++m)
						v2.insert(a->match2[nbr->prey[m]->idx]);
				if(direction == 0 || direction == -1)
					for(m=0;m<nbr->predators.size();++m)
						v2.insert(a->match2[nbr->predators[m]->idx]);
				
				
				//find neighbors shared by species j and k
				v.clear();
				set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(), back_inserter(v));

				v.erase(remove(v.begin(), v.end(), -1), v.end());

				cout << ":";

				totals[0] += v.size();
				totals[1] += v1.size();
				totals[2] += v2.size();
				//print out fraction of interactions shared

				if(v1.size()!=0 && v2.size()!=0){
					cout << v.size() / float(v1.size()); cout << ","; cout << v.size()/float(v2.size());
				}else{
					if(v2.size()!=0){
						cout << 1; cout << ","; cout << v.size()/float(v2.size());
					}else if(v1.size()!=0){
						cout << v.size() / float(v1.size()); cout << ","; cout << 1;
					}else{
						cout << 1; cout << ","; cout << 1;
					}
				}

			}else{
				//Print out alignment of non-aligned nodes
				if (j != -1){

					nbr = n1.nodes[j];
					if(direction == 0 || direction == 1)
						totals[1]+=nbr->prey.size();
					if(direction == 0 || direction == -1)
						totals[1]+=nbr->predators.size();

					cout << n1.roles[j].name;
					cout << ",NULL";
					//Overlap between NULL and species j is set to 0 for both perspectives
					cout << ":0,0";

				}else{

					nbr = n2.nodes[k];
					if(direction == 0 || direction == 1)
						totals[2]+=nbr->prey.size();
					if(direction == 0 || direction == -1)
						totals[2]+=nbr->predators.size();

					cout << "NULL,";
					cout << n2.roles[k].name;
					//Overlap between NULL and species k is set to 0 for both perspectives
					cout << ":0,0";

				}

			}

		if(pairs){
			//Print out pairs distance
			cout << ":";
			cout << distance(a, i);
			cout << ",";
			cout << node_distance(a->matches[i].first, a->matches[i].second, a->dfunc);
		}
		cout <<  ")";
		}

	}
	cout << " ] "; cout << endl;
	cout << "overlap = ("; cout << totals[0] / float(totals[1]);
	cout << ","; cout << totals[0] / float(totals[2]);
	cout << ")"; cout << endl;

}*/



// copy from one alignment to another
void _copy(const void *source, void *dest){
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
    a2->degree = a1->degree;

    a2->fixed_pairs = a1->fixed_pairs;
    a2->set_pairs = a1->set_pairs; 
    a2->unfixed_pairs_A = a1->unfixed_pairs_A;
    a2->unfixed_pairs_B = a1->unfixed_pairs_B;
    a2->doneflag = a1->doneflag;  
}

// copy constructor for an alignment
void * _copy_construct(const void *xp){
	Alignment *a1 = (Alignment *) xp;
	Alignment *a2 = alignment_alloc(a1->match1.size(),a1->match2.size());
	_copy(a1, a2);
	return a2;
}

// destructor for an alignment
void _destroy(void *xp){
	alignment_free((Alignment *) xp);
}

