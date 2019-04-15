#ifndef VMC_ENERGYGRADIENTTARGETFUNCTION_HPP
#define VMC_ENERGYGRADIENTTARGETFUNCTION_HPP

#include "vmc/VMC.hpp"
#include "nfm/NoisyFunction.hpp"

namespace vmc
{

class EnergyGradientTargetFunction: public nfm::NoisyFunctionWithGradient
{
protected:
    VMC &_vmc; // the VMC object containing WaveFunction, Hamiltonian and MCI
    const int _E_Nmc; // number of MC steps for energy calculation
    const int _grad_E_Nmc; // number of MC steps for gradient calculation
    const double _lambda_reg; // vp regularization factor

public:
    EnergyGradientTargetFunction(VMC &vmc, int E_Nmc, int grad_E_Nmc, bool useGradErr, double lambda_reg = 0.):
            nfm::NoisyFunctionWithGradient(vmc.getNVP(), useGradErr), _vmc(vmc),
            _E_Nmc(E_Nmc), _grad_E_Nmc(grad_E_Nmc), _lambda_reg(lambda_reg) {}

    ~EnergyGradientTargetFunction() final = default;

    // NoisyFunctionWithGradient implementation
    nfm::NoisyValue f(const std::vector<double> &vp) final;
    void grad(const std::vector<double> &vp, nfm::NoisyGradient &grad) final;
    nfm::NoisyValue fgrad(const std::vector<double> &vp, nfm::NoisyGradient &grad) final;
};
} // namespace vmc

#endif
