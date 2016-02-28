
// autotools
//#include <config.h>

// c++ includes
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// local includes
#include <common.hpp>
#include <network.hpp>

// namespaces
using namespace std;

Network read_network(char *filename, char separator)
{
  ifstream inFile;
  string line, pred, prey;
  int pred_i, prey_i;
  
  Network N;
  N.name = string(filename);

  inFile.open(filename);
  while (getline(inFile, line)){
    istringstream linestream(line);

    getline(linestream, pred, separator);
    getline(linestream, prey, separator);

    if(N.node_i.count(pred) == 0){
      pred_i = N.nodes.size();
      N.node_i[pred] = pred_i;
      Node * n = new Node;
      n->name = pred;
      N.nodes.push_back(n);
    }else
      pred_i = N.node_i[pred];

    if(N.node_i.count(prey) == 0){
      prey_i = N.nodes.size();
      N.node_i[prey] = prey_i;
      Node * n = new Node;
      n->name = prey;
      N.nodes.push_back(n);
    }else
      prey_i = N.node_i[prey];

    N.nodes[pred_i]->prey.push_back(N.nodes[prey_i]);
    N.nodes[prey_i]->predators.push_back(N.nodes[pred_i]);
  }
  inFile.close();

  return N;
}

void read_networks(char separator, Network& A, Network& B)
{
  string line, pred, prey;
  int pred_i, prey_i;
  
  Network *N = &A;

  A.name = string("stdin");
  B.name = string("stdin");

  while (getline(cin,line)){
    if(line == string("///"))
      N = &B;
    else{
    istringstream linestream(line);

    getline(linestream, pred, separator);
    getline(linestream, prey, separator);

    if(N->node_i.count(pred) == 0){
      pred_i = N->nodes.size();
      N->node_i[pred] = pred_i;
      Node * n = new Node;
      n->name = pred;
      N->nodes.push_back(n);
    }else
      pred_i = N->node_i[pred];

    if(N->node_i.count(prey) == 0){
      prey_i = N->nodes.size();
      N->node_i[prey] = prey_i;
      Node * n = new Node;
      n->name = prey;
      N->nodes.push_back(n);
    }else
      prey_i = N->node_i[prey];

    N->nodes[pred_i]->prey.push_back(N->nodes[prey_i]);
    N->nodes[prey_i]->predators.push_back(N->nodes[pred_i]);
  }
  }
}

void print_network(Network N)
{
  for(unsigned int i=0;i<N.nodes.size();++i)
    for(unsigned int j=0;j<N.nodes[i]->prey.size();++j)
      cout << N.nodes[i]->name << " " << N.nodes[i]->prey[j]->name << endl;
}

// TODO: return a list of the degree-th neighbors to a focal node
set<Node *> neighbors(Network *N, Node *n, unsigned int degree, int direction=0){
    unsigned int i;
    set<Node *> nbrs;
    if(degree == 0){
        nbrs.insert(n);
    }else{
        if(degree == 1){
            if(direction == 0 || direction == 1)
                for(i=0;i<n->prey.size();++i)
                    nbrs.insert(n->prey[i]);
            if(direction == 0 || direction == -1)
                for(i=0;i<n->predators.size();++i)
                    nbrs.insert(n->predators[i]);
        }else{
            cerr << "Sorry, neighborhoods of greater than one degree away are not yet implemented.\n" << endl;
            exit(1);
        }
    }
    return nbrs;
}
