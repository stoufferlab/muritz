
#ifndef __SIMULATED_ANNEALING_HPP
#define __SIMULATED_ANNEALING_HPP

#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include "common.hpp"
#include "anneal.hpp"

double role_euclidean(Role *r1, Role *r2);
double role_correlation(Role *r1, Role *r2);
double role_chisquared(Role *r1, Role *r2);
double role_mahalanobis(Role *r1, Role *r2);

double node_distance(int i, int j);

extern void precompute(unsigned int degree, double (*dfunc) (Role*,Role*));

extern anneal_params_t alignment_params(const gsl_rng *r,
                                        const Alignment *a,
                                        double initialTemperature,
                                        double coolingFactor,
                                        double minTemperature,
                                        int stepsPerTemperature,
                                        int maxUseless,
                                        double acceptanceFraction);

double alignment_get_energy(void *xp);
double alignment_energy_scratch(void *xp);
void   alignment_energy_setup(Alignment *a);
void   alignment_expand_core(void *xp);
void   alignment_rebase_energy(void *xp);
double alignment_propose_step(void *xp, const gsl_rng *r);
void   alignment_commit_step(void *xp);
void   alignment_step(void *xp, const gsl_rng *r);
void   alignment_print(void *xp);
void   alignment_print_pairs(void *xp);
void   overlap_pairs(void *xp, bool pairs, int direction);
void   print_energy(void *xp, cost_func_flag_t cost_function, long degree);

void _copy(const void *source, void *dest);
void _copy_core(const void *source, void *dest);

#endif
