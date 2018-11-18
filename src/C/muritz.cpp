
// autotools
//#include <config.h>

#include <Python.h>
// c++ header files
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>
#include <utility> 
#include <sstream> 

// gsl header files
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_ieee_utils.h>

// for option parsing
#include <unistd.h>

// my header files
#include "common.hpp"
#include "alignment.hpp"
#include "network.hpp"
#include "roles.hpp"
#include "simulated_annealing.hpp"

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

char* muritz(int argc, char *argv[], string net1, string net1_roles, string net2, string net2_roles, string sset_pairs)
{
    // relevant parameters for simulated annealing
    void (*printfunc)(void*) = NULL;
    bool pairs = false, randomstart = false, bipartite=false;
    double iters_fixed_T = 1.0;
    double t_initial = -1.0;
    double mu_t = 1.001;
    double t_min = 1E-7;
    long degree = 0;
    long cost_function = 2;
    int overlap=2;
    // set the above parameters with command line options
    int opt;
    while((opt = getopt(argc, argv, "bvprn:t:c:m:k:l:o:u:")) != -1) {
        switch (opt) {
            case 'b':
                bipartite = true;
                break;
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
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    
    //gsl_rng_set(r, 3716720408280557714);// For testing only.
    
    // read in two files of networks
    n1.bipartite = bipartite;
    n2.bipartite = bipartite;
    read_alignment_data(' ',net1, net1_roles, net2, net2_roles, n1,n2);
    
    //store fixed pairs into a vector
    vector< pair<int, int> > set_pairs;
    vector< int > fixed_pairs;  
    stringstream ss(sset_pairs);
    string pairline;
    int p1_i, p2_i;
    
    //make copies of the species list
    map<string, int> temp_map1; 
    map<string, int> temp_map2; 
    temp_map1.insert(n1.node_i.begin(), n1.node_i.end()); 
    temp_map2.insert(n2.node_i.begin(), n2.node_i.end()); 
    
    while(getline(ss, pairline)) {
        stringstream ssp(pairline);
        string p1, p2;
        ssp >> p1 >> p2;
        map<string, int>::iterator it1 = temp_map1.find(p1);
        map<string, int>::iterator it2 = temp_map2.find(p2);
        if(it1 != temp_map1.end() && it2 != temp_map2.end()) { //if pair is in list
            p1_i = temp_map1[p1];
            p2_i = temp_map2[p2];
            set_pairs.push_back(make_pair(p1_i, p2_i));
            fixed_pairs.push_back(p2_i);
            //delete inserted pair
            temp_map1.erase(it1);
            temp_map2.erase(it2);
        }
    }
    
    // set up the alignment between networks
    Alignment * alignment = setup_alignment(set_pairs);
    
    if(randomstart)
        randomize_alignment(r,alignment);
    
    // decide on what the node-to-node distance function is
    if(cost_function == 0) {
        alignment->dfunc = &role_euclidean_distance;
    } else if(cost_function == 1) {
        alignment->dfunc = &role_correlation;
    } else {
        alignment->dfunc = &role_chisquared;
    }
    
    
    // assign simulated annealing parameters to pass to the function below
    alignment->degree = degree;
    alignment->fixed_pairs = fixed_pairs;
    alignment->set_pairs = set_pairs;
    alignment->doneflag = false;
    
    // Do all the precomputation. Has to be before setting params, as that may make some steps to find the temperature.
    precompute(alignment->degree, alignment->dfunc);
    
    // set up the simulated annealing parameters
    anneal_params_t params = alignment_params(r, alignment, t_initial, mu_t, t_min, iters_fixed_T);
    
    // Set up the initial energy to be altered as the annealing goes on.
    alignment_energy_setup(alignment);
    
    // use simulated annealing to find an optimal alignment
    // print out all of the incremental steps in the the optimization
    anneal(alignment,
           params,
           alignment_energy,
           alignment_propose_step,
           alignment_commit_step,
           printfunc,
           r);
    
    // print out the "optimal" alignment
    alignment->doneflag = true;  
    //alignment_print_json(alignment, true, pairs);
    if(overlap!=0 && overlap!=-1 && overlap!=1){
        
        if(!pairs){
            cout << "optimal ="; alignment_print(alignment);
        }else{
            cout << "optimal ="; alignment_print_pairs(alignment);
        }
        print_energy(alignment, cost_function, degree);
        
    }else{
        //overlap_pairs(alignment, pairs, overlap);
        print_energy(alignment, cost_function, degree);
    }
    
    // free allocated memory
    alignment_free(alignment);
    gsl_rng_free(r);
    
    return "Hello!";
}

/*
char* muritz(char* i) {
    return("Hello!");
}
*/
