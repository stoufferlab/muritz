#include <Python.h>
#include "muritz.h"
#include <iostream>
#include <sstream> 
#include <string>
#include <vector> 

using namespace std; 

static PyObject * muritz_wrapper(PyObject * self, PyObject * args)
{

  char * input; 
  char * cnet1; 
  char * cnet1_roles; 
  char * cnet2; 
  char * cnet2_roles; 
  char * cset_pairs; 
  char * result;
  PyObject * ret;

  // parse arguments
  if (!PyArg_ParseTuple(args, "ssssss", &input, &cnet1, &cnet1_roles, &cnet2, &cnet2_roles, &cset_pairs)) {
    return NULL;
  }
    
  //process input into argv arguments
  vector<char *> argv; 
  int argc = 0; 
  string input_s(input);
  istringstream iss(input_s); 
  string token; 
  while(iss >> token) {
    char *arg = new char[token.size() + 1]; 
    copy(token.begin(), token.end(), arg); 
    arg[token.size()] ='\0';
    argv.push_back(arg); 
    argc++; 
  }
  
  argv.push_back(0); 
  char** argv_p = &argv[0];

  //convert cstrings into strings
  string net1(cnet1); 
  string net1_roles(cnet1_roles); 
  string net2(cnet2);
  string net2_roles(cnet2_roles);
  string set_pairs (cset_pairs); 

  // run the actual function
  result = muritz(argc, argv_p, net1, net1_roles, net2, net2_roles, set_pairs);

  //free space
  for(size_t i = 0; i < argv.size(); i++) {
    delete[] argv[i]; 
  }

  // build the resulting string into a Python object.
  return PyString_FromString(result);
}


static PyMethodDef MuritzMethods[] = {
 { "muritz", muritz_wrapper, METH_VARARGS, "Say hello" },
 { NULL, NULL, 0, NULL }
};


PyMODINIT_FUNC initmuritzex(void)
{
    (void) Py_InitModule("muritzex", MuritzMethods);
}
