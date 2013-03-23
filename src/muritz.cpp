
// autotools
#include <config.h>

// c++ header files
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

// gsl header files
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_ieee_utils.h>

// my header files
#include <common.hpp>
#include <read_roles.hpp>

// namespaces
using namespace std;

int main(int argc, char *argv[])
{
  read_roles(argv[1],' ');
  cout << endl;
  read_roles(argv[2],' ');
  return 0;
}
