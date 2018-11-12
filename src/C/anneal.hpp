#ifndef ANNEALING_H
#define ANNEALING_H

typedef struct {
    double initialTemperature;
    double coolingFactor;
    double minTemperature;
    int stepsPerTemperature;
} anneal_params_t;

typedef double (*anneal_get_energy_t) (void *alignment);
typedef double (*anneal_propose_step_t) (void *alignment);
typedef void (*anneal_make_step_t) (void *alignment);

/*Takes an alignment, parameters for the simulated annealing, and functions to
  interact with the alignment. Performs simualted annealing on the alignment, in place.*/
extern void anneal(void *alignment,
                   anneal_params_t params,
                   anneal_get_energy_t getEnergy,
                   anneal_propose_step_t proposeStep,
                   anneal_make_step_t makeStep);

#endif
