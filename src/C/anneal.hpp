#ifndef ANNEALING_H
#define ANNEALING_H

typedef struct {
    double initialTemperature;
    double coolingFactor;
    double minTemperature;
    int stepsPerTemperature;
    int maxUseless;// The number of temperatures with very low acceptance rate since the last change in the best solution that must be encountered before termination.
    double acceptanceFraction;// The fraction of steps that a temperature must accept fewer than to be considered useless.
} anneal_params_t;

typedef double (*anneal_get_energy_t) (void *alignment);
typedef double (*anneal_propose_step_t) (void *alignment, const gsl_rng *rng);
typedef void (*anneal_make_step_t) (void *alignment);
typedef void (*anneal_print_t) (void *alignment);

/*Takes an alignment, parameters for the simulated annealing, and functions to
  interact with the alignment. Performs simulated annealing on the alignment, in place.*/
extern void anneal(void *alignment,
                   anneal_params_t params,
                   anneal_get_energy_t getEnergy,
                   anneal_propose_step_t proposeStep,
                   anneal_make_step_t makeStep,
                   anneal_print_t printFunc,
                   const gsl_rng *rng);

#endif
