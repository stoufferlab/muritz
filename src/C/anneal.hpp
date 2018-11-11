#ifndef ANNEALING_H
#define ANNEALING_H

struct anneal_params_t {
    double initialTemperature;
    double coolingFactor;
    double minTemperature;
    int stepsPerTemperature;
};

typedef double (*anneal_get_energy_t) (void *alignment);
typedef double (*anneal_propose_step_t) (void *alignment);
typedef void (*anneal_make_step_t) (void *alignment);

extern void anneal(void *alignment,
                   anneal_params_t params,
                   anneal_get_energy_t getEnergy,
                   anneal_propose_step_t proposeStep,
                   anneal_make_step_t makeStep);

#endif
