
#ifndef __ALIGNMENT_HPP
#define __ALIGNMENT_HPP

#include "common.hpp"
#include <vector> 
#include <string>
#include <utility>  

void read_roles(string line, char separator, Network& N);
void read_links(string line, char separator, Network& N);
void read_alignment_data(char separator,string net1, string net1_roles, string net2, string net2_roles, Network& A, Network& B);

Alignment * setup_alignment(vector< pair<int, int> >);
void randomize_alignment(const gsl_rng *r,Alignment *a);

void alignment_print_json(void *xp,bool energy,bool pairs);

Alignment * alignment_alloc(size_t,size_t);
void alignment_free(Alignment *a);

#endif
