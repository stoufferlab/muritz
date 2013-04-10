#ifndef __SIMULATED_ANNEALING_HPP
#define __SIMULATED_ANNEALING_HPP

#include <gsl/gsl_rng.h>
#include <common.hpp>

double role_distance(Role r1, Role r2);
double role_distance(Role r1);
double role_similarity(Role r1, Role r2);
double network_distance(Network net1, Network net2);

double alignment_energy(void *xp);
void   alignment_step(const gsl_rng * r, void *xp, double step_size);
double alignment_distance(void *xp, void *yp);
void   alignment_print(void *xp);

void   alignment_copy(void *source, void *dest);
void * alignment_copy_construct(void *xp);
void   alignment_destroy(void *xp);

#endif