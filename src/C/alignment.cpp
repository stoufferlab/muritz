
// autotools
//#include <config.h>

// c++ header files
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <utility> 

// gsl header files
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

// local includes
#include "common.hpp"
#include "alignment.hpp"
#include "simulated_annealing.hpp"

// namespaces
using namespace std;

// externs
extern Network n1;
extern Network n2;

void read_roles(string input, char separator, Network& N) {
    istringstream in(input);
    string line; 
    bool firstline = true; 
    bool set = false; 
    
    int ncols = 0; 
    string item; 
    char t[1024]; 
    
    while(getline(in, line)) {
        istringstream linestream(line);
        getline(linestream, item, separator);
        Role R;
        R.name = item;
        
        if(firstline){
            firstline = false;
            
            ncols = 1;
            while(getline(linestream, item, separator)){
                Position P;
                sprintf(t,"%i",ncols-1); P.name = string(t);
                P.frequency = strtof(item.c_str(),NULL);
                R.f.push_back(P);      
                ncols++;
            }
        }
        else{
            for(int i=2;i<=ncols;++i){
                getline(linestream, item, separator);
                Position P;
                sprintf(t,"%i",ncols-1); P.name = string(t);
                P.frequency = strtof(item.c_str(),NULL);
                R.f.push_back(P);
            }
        }
        
        N.roles[N.node_i[R.name]] = R;
    }
}

void read_links(string input, char separator, Network& N) {
    string pred, prey; 
    istringstream in(input); 
    string line; 
    int pred_i, prey_i; 
    bool set = false; 
    
    while(getline(in, line)) {
        istringstream linestream(line);
        getline(linestream, pred, separator);
        getline(linestream, prey, separator);
        
        if(N.node_i.count(pred) == 0){
            pred_i = N.nodes.size();
            N.node_i[pred] = pred_i;
            
            Node * n = new Node;
            n->name = pred;
            n->idx = pred_i;
            N.nodes.push_back(n);
            
            Role R;
            N.roles.push_back(R);
            n->side = 0;
        
        }else 
            pred_i = N.node_i[pred];
        
        if(N.node_i.count(prey) == 0){
            prey_i = N.nodes.size();
            N.node_i[prey] = prey_i;
            
            Node * n = new Node;
            n->name = prey;
            n->idx = prey_i;
            N.nodes.push_back(n);
            
            Role R;
            N.roles.push_back(R);
            
            if(N.bipartite){
                n->side = 1;
            } else {
                n->side = 0;
            }
        
        }else
            prey_i = N.node_i[prey];
        
        // add the interactions
        N.nodes[pred_i]->prey.push_back(N.nodes[prey_i]);
        N.nodes[prey_i]->predators.push_back(N.nodes[pred_i]);
    } 
}

void read_alignment_data(char separator,
                         string net1, string net1_roles,
                         string net2, string net2_roles,
                         Network& A, Network& B)
{
    string line, item, pred, prey;
    char t[1024];
    int pred_i, prey_i, ncols = 0;
    bool firstline = true, roles;
    
    vector<string>::iterator it;
    
    A.name = string("Network A");
    B.name = string("Network B");
    
    //network1
    read_links(net1, separator, A); 
    read_roles(net1_roles, separator, A); 
    
    //network2
    read_links(net2, separator, B); 
    read_roles(net2_roles, separator, B); 
}

// setup an alignment structure to manipulate in the SA code
Alignment * setup_alignment(vector< pair<int, int> > set_pairs){
    unsigned int i,j;
    vector<int> unfixed_pairs_A;
    vector<int> unfixed_pairs_B;
    Alignment * a = alignment_alloc(n1.nodes.size(),n2.nodes.size());
    
    // add NULL matches for the nodes in network 1
    for(i=0;i<n1.nodes.size();++i){
        a->matches[i].first = i;
            if (n1.nodes[i]->side == 0){
                unfixed_pairs_A.push_back(i);
            } else {
                unfixed_pairs_B.push_back(i);
            }
    }
    
    // add NULL matches for the nodes in network 2
    for(i=0;i<n2.nodes.size();++i){
        j = i + n1.nodes.size(); // offset based on size of first network
        a->matches[j].second = i;
            if (n2.nodes[i]->side == 0){
                unfixed_pairs_A.push_back(j);
            }else{
                unfixed_pairs_B.push_back(j);
            }
    }
    
    //align the fixed pairs
    int p1_i, p2_i;
    vector<int> fixed_indeces;
    for(i=0; i<set_pairs.size(); i++) {
        p1_i = set_pairs[i].first;
        p2_i = set_pairs[i].second;
        //add fixed matches for the nodes in net1
        a->matches[p1_i].second = p2_i;
        //add fixed matches for the nodes in net2
        a->matches[p2_i + n1.nodes.size()].second =-1;
        //adjust match1 and match2 cheater objects
        a->match1[p1_i] = p2_i;
        a->match2[p2_i] = p1_i;
        //add index to fixed list
        fixed_indeces.push_back(p1_i);
        fixed_indeces.push_back(p2_i + n1.nodes.size());
    }
    
    ptrdiff_t pos;
    //remove fixed pair indices from unfixed pair indices
    sort(fixed_indeces.begin(), fixed_indeces.end());
    for(int j =fixed_indeces.size()-1; j>=0; j--) { //descending
         pos = find(unfixed_pairs_A.begin(), unfixed_pairs_A.end(), fixed_indeces[j]) - unfixed_pairs_A.begin();
         if(pos >= unfixed_pairs_A.size()) {
              pos = find(unfixed_pairs_B.begin(), unfixed_pairs_B.end(), fixed_indeces[j]) - unfixed_pairs_B.begin();
              unfixed_pairs_B.erase(unfixed_pairs_B.begin() + pos);
         } else {
              unfixed_pairs_A.erase(unfixed_pairs_A.begin() + pos);
         }
    }
    
    a->unfixed_pairs_A = unfixed_pairs_A;
    a->unfixed_pairs_B = unfixed_pairs_B;
    
    return a;
}

// randomize an alignment
void randomize_alignment(const gsl_rng *r, Alignment *a){
    unsigned long shuffles = 2*gsl_pow_2(a->unfixed_pairs_A.size()+a->unfixed_pairs_B.size());
    for(unsigned long i=0;i<shuffles;++i) {
        alignment_step(a,r);
    }
}

// print out an alignment in json format
void alignment_print_json(void *xp, bool energy=true, bool pairs=false){
    Alignment * a = (Alignment *) xp;
    unsigned int i;
    int j, k;
    //Role r1, r2;
    
    cout << "{\n";
    cout << "  \"alignment\": {\n";
    
    cout << "    \"energy\": ";
    cout << alignment_get_energy(a) << ",\n";
    
    cout << "    \"pairs\": [\n";
    
    bool first = true;
    for(i=0;i<a->matches.size();++i){
        j = a->matches[i].first;
        k = a->matches[i].second;
        
        // don't print out double NULL matches   
        if(j!=-1 || k!=-1){
            if(!first){
                cout << ",\n";
            }
            
            cout << "      {\n";
            cout << "        \"first\": \"";
            if (j != -1)
                cout << n1.roles[j].name;
            else
                cout << "NULL";
            cout << "\",\n";
            
            cout << "        \"second\": \"";
            if (k != -1)
                cout << n2.roles[k].name;
            else
                cout << "NULL";
            
            if(pairs){
                cout << "\",\n";
                cout << "        \"distance\": ";
                //cout << distance(a, i) << endl; //TODO: Sort this out again.
            }
            
            cout << "      }";
            first = false;
        }
    }
    
    // pairs
    cout << "\n    ]\n";
    
    // alignment
    cout << "  }\n}\n";
        cout << a->matches[i].first << "   " << a->matches[i].second << endl; 
}


// allocate a new alignment
Alignment * alignment_alloc(size_t foo, size_t bar){
    unsigned int i;
    Alignment * a = new Alignment;
    for(i=0;i<foo;++i){
        a->matches.push_back(pair<int,int>(-1,-1));
        a->match1.push_back(-1);
    }
    for(i=0;i<bar;++i){
        a->matches.push_back(pair<int,int>(-1,-1));
        a->match2.push_back(-1);
    }
    return a;
}

// free all memory associated with an alignment
void alignment_free(Alignment * a) {
    delete a;
}
