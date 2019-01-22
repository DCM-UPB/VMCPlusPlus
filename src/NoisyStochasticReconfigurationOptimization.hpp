#ifndef NOISY_STOCHASTIC_RECONFIGURATION_OPTIMIZATION
#define NOISY_STOCHASTIC_RECONFIGURATION_OPTIMIZATION

#include "WFOptimization.hpp"
#include "NoisyStochasticReconfigurationTargetFunction.hpp"
#include "DynamicDescent.hpp"

#include "MCIntegrator.hpp"



class NoisyStochasticReconfigurationOptimization: public WFOptimization{

private:
    long _Nmc;

public:
    NoisyStochasticReconfigurationOptimization(WaveFunction * wf, Hamiltonian * H, const long &Nmc, MCI * mci): WFOptimization(wf, H, mci){
        _Nmc = Nmc;
    }
    virtual ~NoisyStochasticReconfigurationOptimization(){}

    // optimization
    void optimizeWF(){
        // create targetfunction
        NoisyStochasticReconfigurationTargetFunction * targetf = new NoisyStochasticReconfigurationTargetFunction(_wf, _H, getMCI(), _Nmc, 0., true);
        // declare the Dynamic Descent object
        DynamicDescent * ddesc = new DynamicDescent(targetf);
        // allocate an array that will contain the wave function variational parameters
        double * wfpar = new double[_wf->getNVP()];
        // get the variational parameters
        _wf->getVP(wfpar);
        // set the actual variational parameters as starting point for the Conjugate Gradient algorithm
        ddesc->setX(wfpar);
        // find the optimal parameters by minimizing the energy with the Conjugate Gradient algorithm
        ddesc->findMin();
        // set the found parameters in the wave function
        ddesc->getX(wfpar);
        _wf->setVP(wfpar);
        // free memory
        delete[] wfpar;
        delete ddesc;
        delete targetf;
    }
};



#endif
