#ifndef VMC_STOCHASTICRECONFIGURATIONOPTIMIZATION_HPP
#define VMC_STOCHASTICRECONFIGURATIONOPTIMIZATION_HPP

#include "nfm/DynamicDescent.hpp"
#include "vmc/StochasticReconfigurationTargetFunction.hpp"
#include "vmc/WFOptimization.hpp"

#include "mci/MCIntegrator.hpp"

namespace vmc
{

class StochasticReconfigurationOptimization: public WFOptimization
{
private:
    const int _Nmc;
    const double _stepSize;
    const bool _flag_dgrad;

public:
    StochasticReconfigurationOptimization(WaveFunction * wf, Hamiltonian * H, int Nmc, mci::MCI * mci, double stepSize = 1., bool flag_dgrad = false):
            WFOptimization(wf, H, mci), _Nmc(Nmc), _stepSize(stepSize), _flag_dgrad(flag_dgrad) {}

    ~StochasticReconfigurationOptimization() override = default;

    // optimization
    void optimizeWF() override
    {
        // create targetfunction
        StochasticReconfigurationTargetFunction targetf(_wf, _H, getMCI(), _Nmc, 0., _flag_dgrad);
        // declare the Dynamic Descent object
        nfm::DynamicDescent ddesc(targetf.getNDim(), nfm::DDMode::SGDM, false, _stepSize);
        // allocate an array that will contain the wave function variational parameters
        double wfpar[_wf->getNVP()];
        // get the variational parameters
        _wf->getVP(wfpar);
        // set the actual variational parameters as starting point for the Conjugate Gradient algorithm
        ddesc.setX(wfpar);
        // find the optimal parameters by minimizing the energy with the Conjugate Gradient algorithm
        ddesc.findMin(targetf);
        // set the found parameters in the wave function
        _wf->setVP(ddesc.getX().data());
    }
};
} // namespace vmc

#endif
