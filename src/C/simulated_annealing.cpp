
// autotools
//#include <config.h>

// c++ header files
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <iostream>
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
#include "pca.hpp"
#include "simulated_annealing.hpp"

// namespaces
using namespace std;

// global constants
// desired_initial_acceptance is the probability with which a
// maximum-increase-in-energy step is accepted at the start of the annealing.
const double desired_initial_acceptance = 0.5;

// externs
extern Network n1;
extern Network n2;
extern double nullcost;

// globals
gsl_matrix * distance_matrix; // node-to-node distances
double* nulldist1; // distances when unaligned
double* nulldist2; // distances when unaligned

// calculate the role-to-role euclidean distance
double role_euclidean(Role *r1, Role *r2){
	if(r1->name == "NULL" || r2->name == "NULL") {
		return nullcost; // this corresponds to complete lack of correlation
	} else {
		double distance = 0;
		for(unsigned int i=0;i<r1->f.size();++i){
			distance += (r1->f[i].frequency - r2->f[i].frequency ) * (r1->f[i].frequency - r2->f[i].frequency);
		}
		return sqrt(distance);
	}
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

// Calculate the role-to-role distance using Mahalanobis distance.
// This can be calculated by using Principal Coordinate Analysis to transform all the roles,
// normalising them to have unit variance in every dimension,
// and taking the Euclidean distance in the new space.
double role_mahalanobis(Role *r1, Role *r2){
    //The way muritz is currently implemented, only one distance function will ever be called.
    //So it's safe to just go and transform all the role data.
    //If this ever changes in the future, this function will need rewriting.
    //But for now, this is the fastest way to do it.
    static bool are_roles_transformed = false;
    if(!are_roles_transformed) {
        vector<Network*> nets = {&n1, &n2};
        pca_norm_roles(nets);
        are_roles_transformed=true;
    }
    //After transformation, the set of dimensions are linearly uncorrelated.
    //In addition the data has been normalised to have unit variance in every dimension.
    return role_euclidean(r1, r2);
}

// calculate a full role-to-role distance matrix to speed up the SA
static void prepare_distance_matrix(double (*dfunc) (Role*,Role*)){
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
    
    
    /*
    //Print the distance matrix.
    for(unsigned int i = 0; i < n1.nodes.size(); i++) {
        cout << "nulldist1[" << n1.roles[i].name << "] = " << nulldist1[i] << endl;
    }
    for(unsigned int i = 0; i < n2.nodes.size(); i++) {
        cout << "nulldist2[" << n2.roles[i].name << "] = " << nulldist2[i] << endl;
    }
    for(unsigned int i = 0; i < n1.nodes.size(); i++) {
        for(unsigned int j = 0; j < n2.nodes.size(); j++) {
            cout << "dist[" << n1.roles[i].name << "][" << n2.roles[j].name << "] = " << node_distance(i, j) << endl;
        }
    }
    //*/
}

// calculate the nth degree neighbor lists to speed up the neighborhood-based SA
static void prepare_neighbor_data(unsigned int degree){
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
    
    /*
    //Print neighbour data
    cout << "Neighbour data for degree " << degree << ":" << endl;
    cout << "Network 1:" << endl;
    for(unsigned int i = 0; i < n1.nodes.size(); i++) {
        cout << n1.roles[i].name << ":";
        for(nbrs_it = n1.nodes[i]->neighbors[degree].begin(); nbrs_it != n1.nodes[i]->neighbors[degree].end(); nbrs_it++){
            cout << " " << n1.roles[(*nbrs_it)->idx].name;
        }
        cout << endl;
    }
    cout << "Network 2:" << endl;
    for(unsigned int i = 0; i < n2.nodes.size(); i++) {
        cout << n2.roles[i].name << ":";
        for(nbrs_it = n2.nodes[i]->neighbors[degree].begin(); nbrs_it != n2.nodes[i]->neighbors[degree].end(); nbrs_it++){
            cout << " " << n2.roles[(*nbrs_it)->idx].name;
        }
        cout << endl;
    }
    //*/
}

void precompute(unsigned int degree, double (*dfunc) (Role*,Role*)) {
    for(unsigned int d = 0; d <= degree; d++) {
        prepare_neighbor_data(d);
    }
    prepare_distance_matrix(dfunc);
}

// calculate the distance between the roles of two nodes
double node_distance(int i, int j) {
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
    if(delta == 0) return;
    
    map<pair<int, int>, int>::iterator pos = a->matchesContributing.find(make_pair(i, j));
    
    if(pos != a->matchesContributing.end()) {
        // if it already exists in the map
        (*pos).second += delta;// add one
        if((*pos).second == 0) {
            // if there are no copies left in the map
            a->matchesContributing.erase(pos);// remove it
            // this removal ensures that the map size stays linear in the number of nodes in the networks
            // each node can be in the map only twice, once with what it's matched with and once with NULL
        }
    } else {
        // it doesn't yet exist
        a->matchesContributing.insert(make_pair(make_pair(i, j), delta));// insert it
    }
    
    //adjust energy
    a->energy += node_distance(i, j) * delta;
}

// add a contributing match delta to the map, and adjust proposed energy accordingly
static void adjust_proposed_deltas(Alignment *a, int i, int j, int delta) {
    if(delta == 0) return;
    
    //cout << "Adjusting proposed delta between " << (i==-1?"NULL":n1.roles[i].name) << " and " << (j==-1?"NULL":n2.roles[j].name) << " by " << delta << endl;
    
    if(a->proposedContributionDeltas.count(make_pair(i, j))) {
        // if it already exists in the map
        a->proposedContributionDeltas[make_pair(i, j)] += delta;// add delta
    } else {
        // it doesn't yet exist
        a->proposedContributionDeltas.insert(make_pair(make_pair(i, j), delta));// insert it
        // don't both deleting the map entry when it gets to zero: the whole deltas map is wiped regularly
    }
    
    //adjust energy
    a->proposedEnergy += node_distance(i, j) * delta;
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
    
    if(i != -1) nbr_i = n1.nodes[i]->neighbors[degree];
    if(j != -1) nbr_j = n2.nodes[j]->neighbors[degree];
    
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
        a->energy = 0.0;
        // sum the cost function across all paired and unpaired nodes
        for(unsigned int i = 0; i < a->matches.size(); i++){
            a->energy += node_distance(a->matches[i].first, a->matches[i].second);
        }
    } else {// degree != 0
        matches_contributing_setup(a);
    }
}

// expands an alignment created by copy core so it contains energy again
// cannot be further stepped, as it is missing fixed pairs and such
void alignment_expand_core(void *xp) {
    // cast the void parameter as an alignment data type
    Alignment * a = (Alignment *) xp;
    
    // fill match1 and match2
    for(unsigned int i = 0; i < a->matches.size(); i++) {
        int i1 = a->matches[i].first;
        int i2 = a->matches[i].second;
        if(i1 != -1) a->match1[i1] = i2;
        if(i2 != -1) a->match2[i2] = i1;
    }
    
    // set up the energy
    alignment_energy_setup(a);
}

// return the energy/cost function of an alignment
double alignment_get_energy(void *xp)
{
    // cast the void parameter as an alignment data type
    Alignment * a = (Alignment *) xp;
    
    return a->energy;
}

// Rebase the energy of an alignment from the map of contributing pairs, to prevent accumulation of floating point error.
void alignment_rebase_energy(void *xp)
{
    // Cast the void parameter as an alignment data type.
    Alignment * a = (Alignment *) xp;
    
    if(a->degree == 0) {
        alignment_energy_setup(a);// Doing this is no slower.
    } else {// degree != 0
        a->energy = 0.0;
        // Use the map of matches contributing to rebuild the energy.
        for(map<pair<int, int>, int>::iterator it = a->matchesContributing.begin(); it != a->matchesContributing.end(); it++) {
            a->energy += (*it).second * node_distance((*it).first.first, (*it).first.second);
        }
    }
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
    
    //cout << "Propose remove match (" << (i==-1?"NULL":n1.roles[i].name) << "," << (j==-1?"NULL":n1.roles[j].name) << "):" << endl;
    
    if(i != -1) nbr_i = n1.nodes[i]->neighbors[degree];
    if(j != -1) nbr_j = n2.nodes[j]->neighbors[degree];
    
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
        a->proposedContributionWipes.insert(make_pair(i, j));
    }
    // If it's between two non-nulls, wipe all copies of each of them paired with null.
    if(i != -1 && j != -1) {
        a->proposedContributionWipes.insert(make_pair( i, -1));
        a->proposedContributionWipes.insert(make_pair(-1,  j));
    }
}

static void propose_apply_wipes(Alignment *a) {
    
    //cout << "Propose apply wipes:" << endl;
    
    for(set<pair<int, int> >::iterator it = a->proposedContributionWipes.begin(); it != a->proposedContributionWipes.end(); it++) {
        
        int i = (*it).first, j = (*it).second;// The two elements of the match to wipe.
        
        // Find how many already exist, so we know how many to wipe.
        // We don't want to accidentally create a zero, so we can't use []
        map<pair<int, int>, int>::iterator pos = a->matchesContributing.find(make_pair(i, j));
        int numToWipe = 0;
        if(pos != a->matchesContributing.end()) {// if it does exist
            numToWipe = (*pos).second;
        }
        
        // change numToWipe by the delta already entered when adding other matches
        numToWipe += a->proposedContributionDeltas[make_pair(i, j)];
        
        // Adjust the deltas by negative numToWipe, deleting them all
        adjust_proposed_deltas(a, i, j, -numToWipe);
    }
    a->proposedContributionWipes.clear();
}

// add match m, using proposed matches, and adjusting proposedContributionDeltas and proposed energy
static void propose_add_match(Alignment *a, int m) {
    unsigned int degree = a->degree;
    int l;
    int i,j;
    set<Node *> nbr_i, nbr_j;
    set<Node *>::iterator nbr_it;
    
    // How many times we must add *this* match to the network, because adjacent matches have added it.
    int numMatches = 0;
    // How many times we must add the match between this match's network 1 node and null.
    int num1ToNull = 0;
    // How many times we must add the match between this match's network 2 node and null.
    int num2ToNull = 0;
    
    // Note which nodes have been swapped in which networks
    int net1_s1 = a->matches[a->p1].first;
    int net1_s2 = a->matches[a->p2].first;
    int net2_s1 = a->matches[a->p1].second;
    int net2_s2 = a->matches[a->p2].second;
    
    i = a->matches[m].first;
    j = a->matches[m].second;
    // Adjust i and j for proposed swap. Note that swap is always of net2 nodes.
    if(m == a->p1) {
        j = net2_s2;
    } else if(m == a->p2) {
        j = net2_s1;
    }
    
    //cout << "Propose add match (" << (i==-1?"NULL":n1.roles[i].name) << "," << (j==-1?"NULL":n1.roles[j].name) << ")" << endl;
    //cout << "Outgoing:" << endl;
    
    if(i != -1) nbr_i = n1.nodes[i]->neighbors[degree];
    if(j != -1) nbr_j = n2.nodes[j]->neighbors[degree];
    
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
                    // who is i's neighbour aligned to, adjusted for proposed swap?
                    // also note that if adjusting for a swap, we don't add the match itself
                    // because it will be added by the other match, and we don't want to double-count
                    bool ignoreMatch = false;
                    if     (nbr == net1_s1) {l = net2_s2; ignoreMatch = true;}
                    else if(nbr == net1_s2) {l = net2_s1; ignoreMatch = true;}
                    else                    l = a->match1[nbr];
                    
                    // if l is not null and is also one of j's neighbours
                    if(l != -1 && nbr_j.count(n2.nodes[l]) != 0) {
                        adjust_proposed_deltas(a, nbr, l, 1);
                        if(!ignoreMatch) numMatches++;
                    }
                    // l is null or is not one of j's neighbours
                    else
                        adjust_proposed_deltas(a, nbr, -1, 1);
                }
            }else{
                // compute the local alignment for all of j's neighbors
                for(nbr_it=nbr_j.begin(); nbr_it!=nbr_j.end(); ++nbr_it){
                    //what is the index of the neighbour we are currently concerned with?
                    int nbr = (*nbr_it)->idx;
                    // who is j's neighbour aligned to, adjusted for proposed swap?
                    // also note that if adjusting for a swap, we don't add the match itself
                    // because it will be added by the other match, and we don't want to double-count
                    bool ignoreMatch = false;
                    if     (nbr == net2_s1) {l = net1_s2; ignoreMatch = true;}
                    else if(nbr == net2_s2) {l = net1_s1; ignoreMatch = true;}
                    else                    l = a->match2[nbr];
                    
                    // if l is not null and is also one of i's neighbours
                    if(l != -1 && nbr_i.count(n1.nodes[l]) != 0) {
                        adjust_proposed_deltas(a, l, nbr, 1);
                        if(!ignoreMatch) numMatches++;
                    }
                    // l is null or is not one of i's neighbours
                    else
                        adjust_proposed_deltas(a, -1, nbr, 1);
                }
            }
            
            // i and j are not null, so find if we must add (i,NULL) or (NULL,j) as well
            // this is independent of which of i and j has more neighbours, so can't be done in the above loop
            for(nbr_it = nbr_i.begin(); nbr_it != nbr_i.end(); nbr_it++) {
                // what is the index of the neighbour we are currently concerned with?
                int nbr = (*nbr_it)->idx;
                // Who is i's neighbour aligned to?
                // Don't adjust for a proposed swap, just ignore those cases.
                // If the other swapped match would add it, it will add it when we call it.
                // So ignore it here to avoid double-counting.
                if     (nbr == net1_s1) continue;
                else if(nbr == net1_s2) continue;
                else                    l = a->match1[nbr];
                
                // l is null
                if(l == -1) {
                    num1ToNull++;
                }
                // l is not neighbour of j and nbr has more neighbours than l (net1 wins ties)
                else if(nbr_j.count(n2.nodes[l]) == 0 && n1.nodes[nbr]->neighbors[degree].size() >= n2.nodes[l]->neighbors[degree].size()) {
                    num1ToNull++;
                }
            }
            
            for(nbr_it = nbr_j.begin(); nbr_it != nbr_j.end(); nbr_it++) {
                //what is the index of the neighbour we are currently concerned with?
                int nbr = (*nbr_it)->idx;
                // Who is j's neighbour aligned to?
                // Don't adjust for a proposed swap, just ignore those cases.
                // If the other swapped match would add it, it will add it when we call it.
                // So ignore it here to avoid double-counting.
                if     (nbr == net2_s1) continue;
                else if(nbr == net2_s2) continue;
                else                    l = a->match2[nbr];
                
                // l is null
                if(l == -1) {
                    num2ToNull++;
                }
                // l is not neighbour of i and nbr has more neighbours than l (net1 wins ties)
                else if(nbr_i.count(n1.nodes[l]) == 0 && n1.nodes[l]->neighbors[degree].size() < n2.nodes[nbr]->neighbors[degree].size()) {
                    num2ToNull++;
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
                // Don't adjust for a proposed swap, just ignore those cases.
                // If the other swapped match would add it, it will add it when we call it.
                // So ignore it here to avoid double-counting.
                if     (nbr == net1_s1) continue;
                else if(nbr == net1_s2) continue;
                else                    l = a->match1[nbr];
                
                // does nbr have more neighbours than its match (net1 wins ties)?
                // its match being null also qualifies
                if(l == -1 || n1.nodes[nbr]->neighbors[degree].size() >= n2.nodes[l]->neighbors[degree].size()) {
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
                // Don't adjust for a proposed swap, just ignore those cases.
                // If the other swapped match would add it, it will add it when we call it.
                // So ignore it here to avoid double-counting.
                if     (nbr == net2_s1) continue;
                else if(nbr == net2_s2) continue;
                else                    l = a->match2[nbr];
                
                // does nbr have more neighbours than its match (net1 wins ties)?
                // its match being null also qualifies
                if(l == -1 || n1.nodes[l]->neighbors[degree].size() < n2.nodes[nbr]->neighbors[degree].size()) {
                    numMatches ++;
                }
            }
        }
        // else i and j are both null, nothing to add
    }
    
    //cout << "Incoming:" << endl;
    
    // add this match itself the number of times other matches would add it
    adjust_proposed_deltas(a,  i,  j, numMatches);
    adjust_proposed_deltas(a,  i, -1, num1ToNull);
    adjust_proposed_deltas(a, -1,  j, num2ToNull);
}

// calculate the proposed energy from the current energy and the proposed swap
// also updates proposed contributing matches if applicable (if degree != 0)
static void update_proposed_energy(Alignment *a) {
    a->proposedEnergy = a->energy;// Set the current energy to build from.
    a->proposedContributionDeltas.clear();// Clear deltas
    if(a->p1 == a->p2) return;// The proposed step does nothing. TODO: Fix stuff so this can't happen in the first place.
    
    if(a->degree != 0) {
        propose_remove_match(a, a->p1);
        propose_remove_match(a, a->p2);
        propose_apply_wipes(a);
        propose_add_match(a, a->p1);
        propose_add_match(a, a->p2);
    } else {
        a->proposedEnergy -= node_distance(a->matches[a->p1].first, a->matches[a->p1].second);
        a->proposedEnergy -= node_distance(a->matches[a->p2].first, a->matches[a->p2].second);
        a->proposedEnergy += node_distance(a->matches[a->p1].first, a->matches[a->p2].second);
        a->proposedEnergy += node_distance(a->matches[a->p2].first, a->matches[a->p1].second);
    }
}


/* Propose a move in the alignment space and return the energy of the resulting alignment. */
double alignment_propose_step(void *xp, const gsl_rng *r)
{
    // Cast the alignment as an alignment.
    Alignment * a = (Alignment *) xp;
    
    // Find the probability that the switch is of two nodes within the first part (if it's bipartite), as opposed to the second part.
    // If it's not bipartite, this is always 1, as a->unfixed_pairs_B.size() == 0
    static float probA = ((float)a->unfixed_pairs_A.size() / (float)(a->unfixed_pairs_A.size()+a->unfixed_pairs_B.size()));
    
    int p1, p2;
    
    do {
        // Pick the pairs to swap.
        if(gsl_rng_uniform(r) < probA){
            p1 = a->unfixed_pairs_A[gsl_rng_uniform_int(r,a->unfixed_pairs_A.size())];
            p2 = a->unfixed_pairs_A[gsl_rng_uniform_int(r,a->unfixed_pairs_A.size())];
        } else {
            p1 = a->unfixed_pairs_B[gsl_rng_uniform_int(r,a->unfixed_pairs_B.size())];
            p2 = a->unfixed_pairs_B[gsl_rng_uniform_int(r,a->unfixed_pairs_B.size())];
        }
    } while(p1 == p2 || a->matches[p1].first == a->matches[p2].first || a->matches[p1].second == a->matches[p2].second);
    // Do-while ensures that the swap will actually change something; we reroll if it won't.
    // Note that because node indices are unique, the second two conditions only trigger if both relevant nodes are null.
    
    // Put the matches to swap into the alignment.
    a->p1 = p1;
    a->p2 = p2;
    
    // Calculate the energy of the new alignment.
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
    
    apply_proposed_deltas(a);
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
                                 int stepsPerTemperature,
                                 int maxUseless,
                                 double acceptanceFraction) {
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
        ae = alignment_get_energy(b);
        unsigned long shuffles = b->unfixed_pairs_A.size() + b->unfixed_pairs_B.size();
        for(unsigned long i = 0; i < shuffles; i++){
            ae2 = ae;
            alignment_step(b,rng);
            ae = alignment_get_energy(b);
            de = abs(ae - ae2);
            mean_de += de;
            max_de = max(max_de, de);
        }
        mean_de = de/double(shuffles);
        
        // p = e^(-deltaE/T)
        // log(p) = -deltaE/T
        // T = -deltaE/log(p)
        params.initialTemperature = - max_de/log(desired_initial_acceptance);
        
        alignment_free(b);
    }
    
    // damping factor for temperature
    params.coolingFactor = coolingFactor;
    
    // minimum temperature
    params.minTemperature = minTemperature;
    
    // maximum number of temperatures with very few acceptances since finding a new global best
    params.maxUseless = maxUseless;
    
    // fraction of acceptances needed to not count as useless
    params.acceptanceFraction = acceptanceFraction;
    
    return params;
}

// search around up to 'degree' connections away from the aligned nodes and calculate the collective alignment there
// the good old fashioned way: do the whole thing from scratch, without using contributingMatches
static double neighbor_distance_scratch(Alignment *a, unsigned int m) {
    unsigned int degree = a->degree;
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
        // save locally to avoid complications later
        nbr_i = n1.nodes[i]->neighbors[degree];
    }
    
    // j is not null
    if(j != -1){
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
                        d += node_distance((*nbr_it)->idx, l);
                    // l is null or is not one of j's neighbors
                    else
                        d += node_distance((*nbr_it)->idx, -1);
                }
            }else{
                // compute the local alignment for all of j's neighbors
                for(nbr_it=nbr_j.begin(); nbr_it!=nbr_j.end(); ++nbr_it){
                    // who is j's neighbor aligned to?
                    l = a->match2[(*nbr_it)->idx];
                    
                    // if l is not null and is also one of i's neighbors
                    if(l != -1 && nbr_i.count(n1.nodes[l]) != 0)
                        d += node_distance(l, (*nbr_it)->idx);
                    // l is null or is not one of j's neighbors
                    else
                        d += node_distance(-1, (*nbr_it)->idx);
                }
            }
        }
        // j is null
        else{
            // all neighbors of i are treated as unaligned
            for(nbr_it=nbr_i.begin(); nbr_it!=nbr_i.end(); ++nbr_it){
                d += node_distance((*nbr_it)->idx, -1);
            }
        }
    }
    // i is null
    else{
        // j is not null
        if(j != -1){
            // save within a local pointer to avoid complications later
            nbr_j = n2.nodes[j]->neighbors[degree];
            
            // all neighbor nodes are treated as unaligned (i is null)
            for(nbr_it=nbr_j.begin(); nbr_it!=nbr_j.end(); ++nbr_it){
                d += node_distance(-1, (*nbr_it)->idx);
            }
        }
        // j is null
        else{
            d = 0; // for goodness sake...
        }
    }
    
    return d;
}

// calculate the weighted distance between two nodes based on the overall alignment, from scratch
// TODO: set this function up so that we can give it a vector of weights across different "neighborness"
static double distance_scratch(Alignment *a, unsigned int i){
    double d;
    if(a->degree == 0)
        d = node_distance(a->matches[i].first, a->matches[i].second);
    else
        d = neighbor_distance_scratch(a,i);
    return d;
}

// calculate the energy/cost function of an alignment, from scratch
double alignment_energy_scratch(void *xp){
	double E = 0;
	
	// cast the void parameter as an alignment data type
	Alignment * a = (Alignment *) xp;
	
	// sum the cost function across all paired and unpaired nodes
	for(unsigned int i=0;i<a->matches.size();++i){
		E += distance_scratch(a, i);
	}
	return E;
}

// print the energy/cost and the normalized energy/cost function of an alignment
// TODO: Sort out the whole printing deal again.
void print_energy(void *xp, cost_func_flag_t cost_function, long degree){
	double E = 0, Epair_nei, Epair_nod, Enorm_nei=0, Enorm_nod=0;
	int j, k, nei1, nei2, norm=0;
	set<Node *> nbr_i;
	
	// cast the void parameter as an alignment data type
	Alignment * a = (Alignment *) xp;
	
	// sum the cost function across all paired and unpaired nodes
	for(unsigned int i=0;i<a->matches.size();++i){
		Epair_nei = distance_scratch(a, i);
		Epair_nod = node_distance(a->matches[i].first, a->matches[i].second);
		E += Epair_nei;
		if (cost_function == cost_func_flag_t::cor){
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
	
	if (cost_function == cost_func_flag_t::cor && nullcost == 1){
		Enorm_nei=Enorm_nei/(double)norm;
		Enorm_nod=Enorm_nod/(double)norm;
		cout << "Energy = " << E << endl;
		cout << "Normalized nodes energy = " << Enorm_nod << endl;
		if (degree==1){
			cout << "Normalized neighbours energy = " << Enorm_nei << endl;
		}
	}else{
		cout << "Energy = " << E << endl;
	}
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
	cout << " ] " << endl;
	
	/*
	//Print matchesContributing.
	cout << "Matches contributing:" << endl;
	for(map<pair<int, int>, int>::iterator it = a->matchesContributing.begin(); it != a->matchesContributing.end(); it++) {
		cout << "(" << ((*it).first.first==-1?"NULL":n1.roles[(*it).first.first].name) << "," << ((*it).first.second==-1?"NULL":n2.roles[(*it).first.second].name) << ") * " << (*it).second << endl;
	}
	//*/
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
			cout << distance_scratch(a, i);
			cout << ",";
			cout << node_distance(a->matches[i].first, a->matches[i].second);
			cout <<  ")";
		}
	}
	cout << " ] " << endl;
}

//TODO: Probably remove -1 and 1 options.
//They create the size 0 v1 and v2, leading to the counter-intuitive overlap of 1.
void overlap_pairs(void *xp, bool pairs, int direction){
	Alignment * a = (Alignment *) xp;
	unsigned int i,m;
	int j, k;
	vector<int> totals(3);
	set<int> v1;
	set<int> v2;
	vector<int> v;
	Node * node;
	
	
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
				node = n1.nodes[j];
				if(direction == 0 || direction == 1)
					for(m=0;m<node->prey.size();++m)
						v1.insert(node->prey[m]->idx);
				if(direction == 0 || direction == -1)
					for(m=0;m<node->predators.size();++m)
						v1.insert(node->predators[m]->idx);
				
				
				//find alignment partners of neighbors of species k
				v2.clear();
				node = n2.nodes[k];
				if(direction == 0 || direction == 1)
					for(m=0;m<node->prey.size();++m)
						v2.insert(a->match2[node->prey[m]->idx]);
				if(direction == 0 || direction == -1)
					for(m=0;m<node->predators.size();++m)
						v2.insert(a->match2[node->predators[m]->idx]);
				
				
				//find neighbors of j that are aligned to neighbours of k
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
					
					node = n1.nodes[j];
					if(direction == 0 || direction == 1)
						totals[1]+=node->prey.size();
					if(direction == 0 || direction == -1)
						totals[1]+=node->predators.size();
					
					cout << n1.roles[j].name;
					cout << ",NULL";
					//Overlap between NULL and species j is set to 0 for both perspectives
					cout << ":0,0";
					
				}else{
					
					node = n2.nodes[k];
					if(direction == 0 || direction == 1)
						totals[2]+=node->prey.size();
					if(direction == 0 || direction == -1)
						totals[2]+=node->predators.size();
					
					cout << "NULL,";
					cout << n2.roles[k].name;
					//Overlap between NULL and species k is set to 0 for both perspectives
					cout << ":0,0";
					
				}
				
			}
			
			if(pairs){
				//Print out pairs distance
				cout << ":";
				cout << distance_scratch(a, i);
				cout << ",";
				cout << node_distance(a->matches[i].first, a->matches[i].second);
			}
			cout <<  ")";
		}
		
	}
	cout << " ] "; cout << endl;
	cout << "overlap = ("; cout << totals[0] / float(totals[1]);
	cout << ","; cout << totals[0] / float(totals[2]);
	cout << ")"; cout << endl;
}



// copy from one alignment to another
void _copy(const void *source, void *dest) {
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
    
    a2->degree = a1->degree;
    
    a2->fixed_pairs = a1->fixed_pairs;
    a2->set_pairs = a1->set_pairs; 
    a2->unfixed_pairs_A = a1->unfixed_pairs_A;
    a2->unfixed_pairs_B = a1->unfixed_pairs_B;
    a2->doneflag = a1->doneflag;
    
    a2->p1 = a1->p1;
    a2->p2 = a1->p2;
    
    a2->energy = a1->energy;
    a2->proposedEnergy = a1->proposedEnergy;
    
    a2->matchesContributing = a1->matchesContributing;
    a2->proposedContributionDeltas = a1->proposedContributionDeltas;
    a2->proposedContributionWipes = a1->proposedContributionWipes;
}

void _copy_core(const void *source, void *dest) {
    Alignment *a1 = (Alignment *) source, *a2 = (Alignment *) dest;
    
    a2->matches = a1->matches;
    a2->degree = a1->degree;
}
