
#ifndef __ALIGNMENT_HPP
#define __ALIGNMENT_HPP

#include <common.hpp>

void read_alignment_data(char separator,Network& A,Network& B);

Alignment * setup_alignment();
void randomize_alignment(const gsl_rng *r,Alignment *a);

void alignment_print_json(void *xp,bool energy,bool pairs);

Alignment * alignment_alloc(size_t,size_t);
void alignment_free(Alignment *a);

#endif
