
#ifndef __SIMULATED_ANNEALING_HPP
#define __SIMULATED_ANNEALING_HPP

#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <common.hpp>

double role_euclidean_distance(Role *r1, Role *r2);
double role_correlation(Role *r1, Role *r2);
double role_chisquared(Role *r1, Role *r2);

double node_distance(int i, int j, double (*dfunc) (Role*, Role*));
double neighbor_distance(int i, int j, double (*dfunc) (Role*, Role*), unsigned int degree);
double distance(Alignment *a, unsigned int i);

void prepare_distance_matrix(double (*dfunc) (Role*, Role*));
void prepare_neighbor_data(unsigned int degree);

gsl_siman_params_t alignment_params(const gsl_rng * r, void *xp);

double alignment_energy(void *xp);
void   alignment_step(const gsl_rng * r, void *xp, double step_size);
double alignment_distance(void *xp, void *yp);
void   alignment_print(void *xp);
void   alignment_print_pairs(void *xp);
void   overlap_pairs(void *xp, bool pairs, int direction);
void print_energy(void *xp, int cost_function, long degree);


void _copy(void *source, void *dest);
void * _copy_construct(void *xp);
void _destroy(void *xp);

#endif
