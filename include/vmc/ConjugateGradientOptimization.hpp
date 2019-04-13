#ifndef VMC_CONJUGATEGRADIENTOPTIMIZATION_HPP
#define VMC_CONJUGATEGRADIENTOPTIMIZATION_HPP

#include "nfm/ConjGrad.hpp"
#include "vmc/EnergyGradientTargetFunction.hpp"
#include "vmc/WFOptimization.hpp"

#include "mci/MCIntegrator.hpp"

namespace vmc
{

class ConjugateGradientOptimization: public WFOptimization
{
private:
    const int _E_Nmc;
    const int _grad_E_Nmc;

public:
    ConjugateGradientOptimization(WaveFunction * wf, Hamiltonian * H, int E_Nmc, int grad_E_Nmc, mci::MCI * mci):
            WFOptimization(wf, H, mci),
            _E_Nmc(E_Nmc), _grad_E_Nmc(grad_E_Nmc) {}

    ~ConjugateGradientOptimization() override = default;

    // optimization
    void optimizeWF() override
    {
        // create targetfunction
        EnergyGradientTargetFunction targetf(_wf, _H, _E_Nmc, _grad_E_Nmc, getMCI(), true);
        // declare the Conjugate Gradient object
        nfm::ConjGrad cjgrad(targetf.getNDim());
        // allocate an array that will contain the wave function variational parameters
        double wfpar[_wf->getNVP()];
        // get the variational parameters
        _wf->getVP(wfpar);
        // set the actual variational parameters as starting point for the Conjugate Gradient algorithm
        cjgrad.setX(wfpar);
        // find the optimal parameters by minimizing the energy with the Conjugate Gradient algorithm
        cjgrad.findMin(targetf);
        // set the found parameters in the wave function
        _wf->setVP(cjgrad.getX().data());
    }
};
} // namespace vmc

#endif
