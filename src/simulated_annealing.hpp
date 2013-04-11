
#ifndef __SIMULATED_ANNEALING_HPP
#define __SIMULATED_ANNEALING_HPP

#include <gsl/gsl_rng.h>
#include <common.hpp>

double role_distance(Role *r1, Role *r2);
double role_similarity(Role *r1, Role *r2);
//double network_distance(Network net1, Network net2);

Alignment * setup_alignment();

double alignment_energy(void *xp);
void   alignment_step(const gsl_rng * r, void *xp, double step_size);
double alignment_distance(void *xp, void *yp);
void   alignment_print(void *xp);

Alignment * alignment_alloc(size_t n);
void alignment_free(Alignment *a);

void _copy(void *source, void *dest);
void * _copy_construct(void *xp);
void _destroy(void *xp);

#endif
