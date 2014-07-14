
// autotools
//#include <config.h>

// c++ header files
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

// gsl header files
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

// local includes
#include <common.hpp>
#include <alignment.hpp>
#include <simulated_annealing.hpp>

// namespaces
using namespace std;

// externs
extern Network n1;
extern Network n2;

void read_alignment_data(char separator, Network& A, Network& B)
{
  string line, item, pred, prey;
  char t[1024];
  int pred_i, prey_i, ncols = 0;
  bool firstline = true, roles;
  
  Network *N = &A;

  A.name = string("Network A");
  B.name = string("Network B");

  roles = false;
  while (getline(cin,line)){
    if(line == string("///")){
      N = &B;
      roles=false;
    }
    else
      if(line == string("###")){
        roles=true;
        firstline = true;
      }
      else
        if(!roles){
          istringstream linestream(line);
          getline(linestream, pred, separator);
          getline(linestream, prey, separator);

          if(N->node_i.count(pred) == 0){
            pred_i = N->nodes.size();
            N->node_i[pred] = pred_i;
            
            Node * n = new Node;
            n->name = pred;
            n->idx = pred_i;
            N->nodes.push_back(n);
            
            Role R;
            N->roles.push_back(R);
          }else
            pred_i = N->node_i[pred];

          if(N->node_i.count(prey) == 0){
            prey_i = N->nodes.size();
            N->node_i[prey] = prey_i;
            
            Node * n = new Node;
            n->name = prey;
            n->idx = prey_i;
            N->nodes.push_back(n);

            Role R;
            N->roles.push_back(R);
          }else
            prey_i = N->node_i[prey];

          // add the interactions
          N->nodes[pred_i]->prey.push_back(N->nodes[prey_i]);
          N->nodes[prey_i]->predators.push_back(N->nodes[pred_i]);
        }
        else{
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
              P.frequency = strtol(item.c_str(),NULL,10);

              R.f.push_back(P);
      
              ncols++;
            }
          }
          else{
            for(int i=2;i<=ncols;++i){
              getline(linestream, item, separator);
              Position P;
              sprintf(t,"%i",ncols-1); P.name = string(t);
              P.frequency = strtol(item.c_str(),NULL,10);

              R.f.push_back(P);
            }
          }

          N->roles[N->node_i[R.name]] = R;
        }
  }
}

// setup an alignment structure to manipulate in the SA code
Alignment * setup_alignment(){
	unsigned int i,j;
	Alignment * a = alignment_alloc(n1.nodes.size(),n2.nodes.size());

 	// add NULL matches for the nodes in network 1
 	for(i=0;i<n1.nodes.size();++i){
 		a->matches[i].first = i;
  }

 	// add NULL matches for the nodes in network 2
 	for(i=0;i<n2.nodes.size();++i){
 		j = i + n1.nodes.size(); // offset based on size of first network
 		a->matches[j].second = i;
  }

 	return a;
}

// randomize an alignment
void randomize_alignment(const gsl_rng *r, Alignment *a){
  unsigned long shuffles = 2*gsl_pow_2(a->matches.size());
  for(unsigned long i=0;i<shuffles;++i)
    alignment_step(r,a,0);
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
  cout << alignment_energy(a) << ",\n";

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
        cout << distance(a, i) << endl;;
      }
      
      cout << "      }";
      first = false;
    }
  }

  // pairs
  cout << "\n    ]\n";

  // alignment
  cout << "  }\n}\n";
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
void alignment_free(Alignment * a){
	a->matches.clear();
  a->match1.clear();
  a->match2.clear();
	delete a;
}
