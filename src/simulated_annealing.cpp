
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

double role_distance(Role r1, Role r2){
	//cout << "role distance (double)...";
	double distance = 0;
	for(int i=0;i<r1.f.size();++i){
		distance += (r1.f[i].frequency - r2.f[i].frequency ) * (r1.f[i].frequency - r2.f[i].frequency);
	}
	//cout << distance << "...done\n";
  	return distance;
}

double role_distance(Role r1){
	//cout << "role distance (single)...";	
	double distance = 0;
	for(int i=0;i<r1.f.size();++i){
		distance += (r1.f[i].frequency) * (r1.f[i].frequency);
	}
	//cout << distance << "...done\n";
  	return distance;
}

double role_similarity(Role r1, Role r2){
  return 1 - role_distance(r1,r2);
}

double alignment_energy(void *xp){
	//cout << "alignment energy...";
	Alignment *a = (Alignment *) xp;
	double E = 0;
	unsigned int i,j,k;

	for(i=0;i<a->matches.size();++i){
		j = a->matches[i].first;
		k = a->matches[i].second;

		if(j != -1){
			if(k != -1){
				E += role_distance(a->net1->roles[j], a->net2->roles[k]);
			}else{
				E += role_distance(a->net1->roles[j]);
			}
		}else{
			if(k != -1){
				E += role_distance(a->net2->roles[k]);
			}
		}
	}

	cout << "energy = " << E << ": ";
	alignment_print(xp);
	cout << endl;
	return E;
}

/* make a move in the alignment space */
void alignment_step(const gsl_rng * r, void *xp, double step_size){
	Alignment *a = (Alignment *) xp;
	int p1, p2, tmp;

	step_size = 0 ; // prevent warnings about unused parameter

	// pick the pairs to swap
	p1 = gsl_rng_uniform_int(r,a->matches.size());
	p2 = gsl_rng_uniform_int(r,a->matches.size());

	// swap the indices for net2
	tmp = a->matches[p1].second;
	a->matches[p1].second = a->matches[p2].second;
	a->matches[p2].second = tmp;

	if(a->matches[p1].first != -1)
		a->match1[a->matches[p1].first] = a->matches[p1].second;

	if(a->matches[p2].first != -1)
		a->match1[a->matches[p2].first] = a->matches[p2].second;

}

double alignment_distance(void *xp, void *yp){
	//cout << "alignment distance...";
	Alignment *a1 = (Alignment *) xp, *a2 = (Alignment *) yp;
	double distance = 0;
  	unsigned int i;

  	for (i=0; i<a1->match1.size();++i){
		distance += ((a1->match1[i] == a2->match1[i]) ? 0 : 1);
  	} 

	//cout << "done\n";
	return distance;
}

void alignment_print(void *xp){
	Alignment *a = (Alignment *) xp;
	unsigned int i,j,k;

	cout << "  [";
	for(i=0;i<a->matches.size();++i){
		j = a->matches[i].first;
		k = a->matches[i].second;
		
		cout << " (";
		if (j != -1)
			cout << a->net1->roles[j].name;
		else
			cout << "NULL";
		cout << ",";
		if (k != -1)
			cout << a->net2->roles[k].name;
		else
			cout << "NULL";
		cout <<  ")";
	}
	cout << " ]  ";
}

void alignment_copy(void *source, void *dest){
	//cout << "copy...";
	
	Alignment *a1 = (Alignment *) source;
	Alignment *a2 = (Alignment *) dest;
	
	a2->net1 = a1->net1;
	a2->net2 = a1->net2;

	a2->matches.resize(a1->matches.size());
	for(unsigned int i=0;i<a1->matches.size();++i)
		a2->matches[i] = a1->matches[i];

	a2->match1.resize(a1->match1.size());
	for(unsigned int i=0;i<a1->match1.size();++i)
		a2->match1[i] = a1->match1[i];
	
	//cout << "done\n";
}

void * alignment_copy_construct(void *xp){
	//cout << "copy construct...";
	Alignment *a1 = (Alignment *) xp;
	Alignment *a2;

	a2->net1 = a1->net1;
	a2->net2 = a1->net2;

	a2->matches = a1->matches;
	a2->match1 = a1->match1;

	cout << "done\n";
	return a2;
}

void alignment_destroy(void *xp){
	//cout << "destroy...";
	Alignment *a = (Alignment *) xp;

	a->matches.clear();
	a->match1.clear();
	cout << "done\n";
}

/*double network_distance(Network net1, Network net2){
  return 3.14159;
}*/