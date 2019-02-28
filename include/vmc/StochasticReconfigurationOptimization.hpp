#ifndef VMC_STOCHASTICRECONFIGURATIONOPTIMIZATION_HPP
#define VMC_STOCHASTICRECONFIGURATIONOPTIMIZATION_HPP

#include "nfm/DynamicDescent.hpp"
#include "vmc/StochasticReconfigurationTargetFunction.hpp"
#include "vmc/WFOptimization.hpp"

#include "mci/MCIntegrator.hpp"



class StochasticReconfigurationOptimization: public WFOptimization
{
private:
    const int _Nmc;
    const double _stepSize;
    const bool _flag_dgrad;

public:
    StochasticReconfigurationOptimization(WaveFunction * wf, Hamiltonian * H, const int &Nmc, MCI * mci, const double stepSize = 1., const bool flag_dgrad = false): WFOptimization(wf, H, mci), _Nmc(Nmc), _stepSize(stepSize), _flag_dgrad(flag_dgrad) {}

    ~StochasticReconfigurationOptimization() override= default;

    // optimization
    void optimizeWF() override{
        // create targetfunction
        StochasticReconfigurationTargetFunction * targetf = new StochasticReconfigurationTargetFunction(_wf, _H, getMCI(), _Nmc, 0., _flag_dgrad);
        // declare the Dynamic Descent object
        DynamicDescent ddesc(targetf, _stepSize);
        // allocate an array that will contain the wave function variational parameters
        double wfpar[_wf->getNVP()];
        // get the variational parameters
        _wf->getVP(wfpar);
        // set the actual variational parameters as starting point for the Conjugate Gradient algorithm
        ddesc.setX(wfpar);
        // find the optimal parameters by minimizing the energy with the Conjugate Gradient algorithm
        ddesc.findMin();
        // set the found parameters in the wave function
        ddesc.getX(wfpar);
        _wf->setVP(wfpar);

        delete targetf;
    }
};



#endif
