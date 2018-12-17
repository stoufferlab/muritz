#include <math.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <limits>
#include <iostream>
#include "anneal.hpp"
//#include "simulated_annealing.hpp"// Only needed if we want to print the alignment energy calculated from scratch.

using namespace std;

// Taking 1 million to be a conservative upper bound on the energy, a difference of one Unit in the Last Place is 1.16415e-10
// EPSILON of 1e-6 should allow plenty of ULPs of error. TODO: Test this nonetheless. Also ensure it's not too large.
#define ENERGY_EPSILON 1e-6

// Compares energy values for equality.
static inline bool energyEqual(double a, double b) {
    return abs(a-b) < ENERGY_EPSILON;
}

/*Return a number in the range [0,1] that is the probability of taking the described step.*/
static double getStepProbability(double oldEnergy, double newEnergy, double temperature)
{
    if(newEnergy < oldEnergy) {
        return 1.0;//If this step improves the cost, always take it.
    } else if(temperature == 0) {
        return 0.0;//Never take steps that don't improve the cost at temperature zero.
    } else {
        return exp((oldEnergy-newEnergy)/temperature);
    }
}



void anneal(void *alignment,
            void *bestAlignment,
            anneal_params_t params,
            anneal_get_energy_t getEnergy,
            anneal_propose_step_t proposeStep,
            anneal_commit_step_t commitStep,
            anneal_rebase_energy_t rebaseEnergy,
            anneal_print_t printFunc,
            anneal_copy_core_t copyCore,
            const gsl_rng *rng)
{
    int numAcceptsNeeded =
        (int) ceil(params.stepsPerTemperature * params.acceptanceFraction);
    
    double temperature = params.initialTemperature;
    double prevTemperature = temperature; //For printing the final temperature.
    double currentEnergy = getEnergy(alignment);
    double nextEnergy;
    double bestEnergy = numeric_limits<double>::infinity();
    int numUseless = 0;
    
    //printf("Steps per temperature: %d\n", params.stepsPerTemperature);
    
    while(temperature != 0 && (numUseless < params.maxUseless || params.maxUseless == 0)) {
        if(temperature < params.minTemperature) temperature = 0;//Do one final run of pure hill-climbing.
        
        //printf("Temperature = %.12lf\n", temperature);
        
        // rebase the energy to prevent accumulation of floating point error
        rebaseEnergy(alignment);
        currentEnergy = getEnergy(alignment);
        
        int numAccepts = 0;
        
        for(int i = 0; i < params.stepsPerTemperature; i++) {
            nextEnergy = proposeStep(alignment, rng);//Propose a step and get its cost.
            
            //Calculate the probability of taking the proposed step.
            double probability = getStepProbability(currentEnergy, nextEnergy, temperature);
            
            //if(printFunc) printf("Current energy = %.12lf, proposed energy = %.12lf", currentEnergy, nextEnergy);
            //if(printFunc) printf(", current energy from scratch = %.12lf", alignment_energy_scratch(alignment));
            
            if(gsl_rng_uniform(rng) < probability) {// Take the step.
                commitStep(alignment);
                
                if(!energyEqual(currentEnergy, nextEnergy)) {
                    //Don't count it if the energies are equal.
                    //Otherwise it will likely flip-flop between
                    //two equal-energy alignments and never terminate early.
                    numAccepts++;
                }
                
                currentEnergy = nextEnergy;
                
                /*if(printFunc) {
                    printf(", taking step. New alignment:\n");
                    printFunc(alignment);
                }*/
            }/* else {// Don't take the step.
                if(printFunc) {
                    printf(", rejecting step.\n");
                }
            }*/
            
            if(currentEnergy < bestEnergy && !energyEqual(currentEnergy, bestEnergy)) {
                copyCore(alignment, bestAlignment);
                bestEnergy = currentEnergy;
                numUseless = 0;
            }
        }
        
        if(printFunc) {
            printf("temperature: %12g, energy: %12g, best energy: %12g, alignment:\n",
                temperature, currentEnergy, bestEnergy);
            printFunc(alignment);
        }
        
        if(numAccepts <= numAcceptsNeeded) numUseless++;
        
        prevTemperature = temperature;
        temperature /= params.coolingFactor;
    }
    
    cout << "final temperature: " << prevTemperature << endl;
}
