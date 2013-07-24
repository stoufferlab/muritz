
#ifndef __SIMULATED_ANNEALING_HPP
#define __SIMULATED_ANNEALING_HPP

#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <common.hpp>

double role_euclidean_distance(Role *r1, Role *r2);
double role_correlation(Role *r1, Role *r2);
void role_chisquared(Role *r1, Role *r2, double& chisq, int& df);
double role_distance(Role *r1, Role *r2);

void prepare_distance_matrix(void);
gsl_siman_params_t alignment_params(void *xp);

double alignment_energy(void *xp);
void   alignment_step(const gsl_rng * r, void *xp, double step_size);
double alignment_distance(void *xp, void *yp);
void   alignment_print(void *xp);

void _copy(void *source, void *dest);
void * _copy_construct(void *xp);
void _destroy(void *xp);

#endif
