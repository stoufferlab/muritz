
// autotools
#include <config.h>

// c++ header files
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// gsl header files
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

// local includes
#include <common.hpp>
#include <simulated_annealing.hpp>

// namespaces
using namespace std;

// externs
extern Network n1;
extern Network n2;

double role_distance(Role *r1, Role *r2){
	double distance = 0;
	for(int i=0;i<r1->f.size();++i){
		if(r2->name != "NULL")
			distance += (r1->f[i].frequency - r2->f[i].frequency ) * (r1->f[i].frequency - r2->f[i].frequency);
		else
			distance += (r1->f[i].frequency) * (r1->f[i].frequency);
	}
  	return distance;
}

double role_similarity(Role *r1, Role *r2){
  return 1 - role_distance(r1,r2);
}

// calculate the energy/cost function of an alignment
double alignment_energy(void *xp){
	//cout << "energying\n";
	Alignment * a = (Alignment *) xp;
	unsigned int i, j, k;
	double E = 0;
	Role r1, r2, null;
	null.name = "NULL";

	for(i=0;i<a->matches.size();++i){
		j = a->matches[i].first;
		k = a->matches[i].second;

		if(j != -1){
			r1 = n1.roles[j];
			if(k != -1){
				r2 = n2.roles[k];
				E += role_distance(&r1,&r2);
			}else{
				E += role_distance(&r1,&null);
			}
		}else{
			if(k != -1){
				r2 = n2.roles[k];
				E += role_distance(&r2,&null);
			}
		}
	}

	//cout << "energy = " << E << ": ";
	//alignment_print(xp);
	//cout << endl;
	return E;
}

/* make a move in the alignment space */
void alignment_step(const gsl_rng * r, void *xp, double step_size){
	//cout << "stepping\n";
	step_size = 0 ; // prevent warnings about unused parameter

	Alignment * a = (Alignment *) xp;

	// pick the pairs to swap
	unsigned int p1 = gsl_rng_uniform_int(r,a->matches.size());
	unsigned int p2 = gsl_rng_uniform_int(r,a->matches.size());

	// swap the indices for net2
	unsigned int tmp = a->matches[p1].second;
	a->matches[p1].second = a->matches[p2].second;
	a->matches[p2].second = tmp;
}

// calculate the distance between two alignments
double alignment_distance(void *xp, void *yp){
	//cout << "distancing\n";
	Alignment *a1 = (Alignment *) xp, *a2 = (Alignment *) yp;
	double distance = 0;
  	for(unsigned int i=0; i<a1->matches.size();++i){
  		// check if each pairwise match is the same
		distance += ((a1->matches[i] == a2->matches[i]) ? 0 : 1);
  	} 
	return distance;
}

// print out an alignment
void alignment_print(void *xp){
	//cout << "printing\n";
	Alignment * a = (Alignment *) xp;
	unsigned int i, j, k;
	Role r1, r2;

	cout << "  [";
	for(i=0;i<a->matches.size();++i){
		j = a->matches[i].first;
		k = a->matches[i].second;

		// don't print out NULL matches		
		if(j!=-1 || k!=-1){
			cout << " (";
			if (j != -1)
				cout << n1.roles[j].name;
			else
				cout << "NULL";
			cout << ",";
			if (k != -1)
				cout << n2.roles[k].name;
			else
				cout << "NULL";
			cout <<  ")";
		}
	}
	cout << " ]  ";
}

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

// copy from one alignment to another
void _copy(void *source, void *dest){
	Alignment *a1 = (Alignment *) source, *a2 = (Alignment *) dest;
	for(unsigned int i=0;i<a1->matches.size();++i){
		a2->matches[i] = a1->matches[i];
	}
}

// copy constructor for an alignment
void * _copy_construct(void *xp){
	Alignment *a1 = (Alignment *) xp;
	Alignment *a2 = alignment_alloc(a1->matches.size());
	_copy(a1, a2);
	return a2;
}

// destructor for an alignment
void _destroy(void *xp){
	alignment_free((Alignment *) xp);
}