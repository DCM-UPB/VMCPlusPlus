#ifndef ADAM_STOCHASTIC_RECONFIGURATION_OPTIMIZATION
#define ADAM_STOCHASTIC_RECONFIGURATION_OPTIMIZATION

#include "WFOptimization.hpp"
#include "NoisyStochasticReconfigurationTargetFunction.hpp"
#include "Adam.hpp"

#include "MCIntegrator.hpp"



class AdamStochasticReconfigurationOptimization: public WFOptimization{

private:
    long _Nmc;
    long _grad_E_Nmc;
    double _lambda;
    double _alpha;
    double _beta1, _beta2;
    double _epsilon;

public:
    AdamStochasticReconfigurationOptimization(WaveFunction * wf, Hamiltonian * H, MCI * mci, const long &Nmc, const long &grad_E_Nmc, const double &lambda = 0., const double &alpha = 0.001, const double &beta1 = 0.9, const double &beta2 = 0.999, const double &epsilon = 10e-8)
        : WFOptimization(wf, H, mci){
        _Nmc = Nmc;
        _grad_E_Nmc = grad_E_Nmc;
        _lambda = lambda;
        _alpha = alpha;
        _beta1 = beta1;
        _beta2 = beta2;
        _epsilon = epsilon;
    }
    virtual ~AdamStochasticReconfigurationOptimization(){}

    // optimization
    void optimizeWF(){
        // create targetfunction
        NoisyStochasticReconfigurationTargetFunction * targetf = new NoisyStochasticReconfigurationTargetFunction(_wf, _H, _Nmc, getMCI(), false);
        // declare the Adam object
        Adam * adam = new Adam(targetf, _alpha, _beta1, _beta2, _epsilon);
        // allocate an array that will contain the wave function variational parameters
        double wfpar[_wf->getNVP()];
        // get the variational parameters
        _wf->getVP(wfpar);
        // set the actual variational parameters as starting point for optimization
        adam->setX(wfpar);
        // find the optimal parameters by minimizing the energy with the Adam algorithm
        adam->findMin();
        // set the found parameters in the wave function
        adam->getX(wfpar);
        _wf->setVP(wfpar);
        // free memory
        delete adam;
        delete targetf;
    }

};



#endif
