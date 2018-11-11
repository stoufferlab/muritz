#include <math.h>
#include <random>
#include "anneal.hpp"

/*Return a number in the range [0,1] that is the probability of taking the described step.*/
static double getStepProbability(double oldEnergy, double newEnergy, double temperature)
{
    if(newEnergy < oldEnergy) {
        return 1.0;
    } else if(temperature == 0) {
        return 0.0;
    } else {
        return exp((oldEnergy-newEnergy)/temperature);
    }
}



extern void anneal(void *alignment,
                   anneal_params_t params,
                   anneal_get_energy_t getEnergy,
                   anneal_propose_step_t proposeStep,
                   anneal_make_step_t makeStep)
{
    double temperature = params.initialTemperature;
    double currentEnergy = getEnergy(alignment);
    double nextEnergy;
    
    std::random_device rd;//To seed the LCG.
    std::minstd_rand gen(rd());//Seed the LCG from rd.
    std::uniform_real_distribution<> runif(0.0,1.0);
    //runif(gen) will generate a random number in the range [0,1)
    
    while(temperature != 0) {
        if(temperature < params.minTemperature) temperature = 0;
        
        for(int i = 0; i < params.stepsPerTemperature; i++) {
            nextEnergy = proposeStep(alignment);
            double probability = getStepProbability(currentEnergy, nextEnergy, temperature);
            
            if(runif(gen) < probability) {
                makeStep(alignment);
            }
        }
    }
}
