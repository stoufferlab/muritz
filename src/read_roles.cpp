
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

// namespaces
using namespace std;

Network read_roles(char *filename, char separator)
{
  ifstream inFile;
  string line,item;
  char t[1024];
  Network N;
  N.name = string(filename);
  bool firstline = true;
  int ncols;

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

    N.roles.push_back(R);
  }
  inFile.close();

  return N;
}
