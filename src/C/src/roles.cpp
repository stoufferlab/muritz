
// autotools
#include <config.h>

// c++ includes
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// local includes
#include <common.hpp>
#include <roles.hpp>

// namespaces
using namespace std;

vector<Role> read_roles(char *filename, char separator, Network N)
{
  ifstream inFile;
  string line,item;
  char t[1024];
  bool firstline = true;
  int ncols=0;

  // start with an empty role vector
  vector<Role> roles(N.nodes.size());

  inFile.open(filename);
  while (getline(inFile, line)){
    istringstream linestream(line);
    Role R;

    if(firstline){
      firstline = false;
      
      getline(linestream, item, separator);
      R.name = item;
      
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
      getline(linestream, item, separator);
      R.name = item;

      for(int i=2;i<=ncols;++i){
        getline(linestream, item, separator);
        Position P;
        sprintf(t,"%i",ncols-1); P.name = string(t);
        P.frequency = strtol(item.c_str(),NULL,10);

        R.f.push_back(P);
      }
    }

    roles[N.node_i[R.name]] = R;
  }
  inFile.close();

  return roles;
}

vector<Role> read_roles(char separator, Network N)
{
  string line,item;
  char t[1024];
  bool firstline = true;
  int ncols=0;

  // start with an empty role vector
  vector<Role> roles(N.nodes.size());

  while (getline(cin, line)){
    istringstream linestream(line);
    Role R;

    if(firstline){
      firstline = false;
      
      getline(linestream, item, separator);
      R.name = item;
      
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
      getline(linestream, item, separator);
      R.name = item;

      for(int i=2;i<=ncols;++i){
        getline(linestream, item, separator);
        Position P;
        sprintf(t,"%i",ncols-1); P.name = string(t);
        P.frequency = strtol(item.c_str(),NULL,10);

        R.f.push_back(P);
      }
    }

    roles[N.node_i[R.name]] = R;
  }

  return roles;
}

void print_roles(Network N){
  unsigned int i,j;
  for(i=0;i<N.nodes.size();++i){
    cout << N.nodes[i]->name << " ";
    for(j=0;j<N.roles[i].f.size()-1;++j)
      cout << N.roles[i].f[j].frequency << " ";
    cout << N.roles[i].f[j].frequency << endl;
  }
}