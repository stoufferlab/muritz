
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

double role_distance(Role r1, Role r2){
  return 0.5;
}

double role_similarity(Role r1, Role r2){
  return 1 - role_distance(r1,r2);
}

double network_distance(Role r1, Role r2){
  return 3.14159;
}