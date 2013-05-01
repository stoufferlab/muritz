
// autotools
#include <config.h>

// c++ header files
#include <iostream>
#include <vector>

// local includes
#include <common.hpp>
#include <alignment.hpp>

// namespaces
using namespace std;

// externs
extern Network n1;
extern Network n2;

// setup an alignment structure to manipulate in the SA code
Alignment * setup_alignment(){
	unsigned int i,j,k;
	Alignment * a = alignment_alloc(n1.roles.size()+n2.roles.size());

  	// add NULL matches for the nodes in network 1
  	for(i=0;i<n1.roles.size();++i){
  		a->matches[i].first = i;
  		a->matches[i].second = -1;
  	}

  	// add NULL matches for the nodes in network 2
  	for(i=0;i<n2.roles.size();++i){
  		j = n1.roles.size() + i;
  		a->matches[j].first = -1;
  		a->matches[j].second = i;
  	}

  	return a;
}

// allocate a new alignment
Alignment * alignment_alloc(size_t n){
	Alignment * a = new Alignment;
	for(unsigned int i=0;i<n;++i)
  		a->matches.push_back(pair<int,int>(-1,-1));
	return a;
}

// free all memory associated with an alignment
void alignment_free(Alignment * a){
	a->matches.clear();
	delete a;
}
