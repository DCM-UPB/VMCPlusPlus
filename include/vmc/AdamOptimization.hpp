#ifndef ADAM_OPTIMIZATION
#define ADAM_OPTIMIZATION

#include "vmc/WFOptimization.hpp"
#include "vmc/StochasticReconfigurationTargetFunction.hpp"
#include "vmc/EnergyGradientTargetFunction.hpp"
#include "nfm/Adam.hpp"

#include "mci/MCIntegrator.hpp"



class AdamOptimization: public WFOptimization{

private:
    long _Nmc;
    bool _useSR;
    bool _useGradientError;
    bool _useAveraging;
    size_t _max_n_const_values;
    double _lambda;
    double _alpha;
    double _beta1, _beta2;
    double _epsilon;

public:
    AdamOptimization(WaveFunction * wf, Hamiltonian * H, MCI * mci, const long &Nmc, const bool useSR = false, const bool useGradientError = false, const size_t &max_n_const_values = 20,
        const bool useAveraging = false, const double &lambda = 0., const double &alpha = 0.001, const double &beta1 = 0.9, const double &beta2 = 0.999, const double &epsilon = 10e-8)
        : WFOptimization(wf, H, mci)
    {
        _Nmc = Nmc;
        _useSR = useSR;
        _useGradientError = useGradientError;
        _max_n_const_values = max_n_const_values;
        _useAveraging = useAveraging;
        _lambda = lambda;
        _alpha = alpha;
        _beta1 = beta1;
        _beta2 = beta2;
        _epsilon = epsilon;
    }
    virtual ~AdamOptimization(){}

    // optimization
    void optimizeWF(){
        // create targetfunction
        NoisyFunctionWithGradient * targetf;
        if (_useSR) targetf = new StochasticReconfigurationTargetFunction(_wf, _H, getMCI(), _Nmc, _lambda, false);
        else targetf = new EnergyGradientTargetFunction(_wf, _H, _Nmc, _Nmc, getMCI(), _lambda);
        // declare the Adam object
        Adam * adam = new Adam(targetf, _useGradientError, _max_n_const_values, _useAveraging, _alpha, _beta1, _beta2, _epsilon);
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
