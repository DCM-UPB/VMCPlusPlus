#ifndef VMC_STOCHASTICRECONFIGURATIONTARGETFUNCTION_HPP
#define VMC_STOCHASTICRECONFIGURATIONTARGETFUNCTION_HPP

#include "vmc/VMC.hpp"
#include "nfm/NoisyFunction.hpp"

namespace vmc
{

class StochasticReconfigurationTargetFunction: public nfm::NoisyFunctionWithGradient
{
protected:
    VMC &_vmc; // the VMC object containing WaveFunction, Hamiltonian and MCI
    const int64_t _E_Nmc; // number of MC steps for energy calculation
    const int64_t _grad_E_Nmc; // number of MC steps for gradient calculation
    const double _lambda_reg; // vp regularization factor

    void _integrate(const double * vp, double * obs, double * dobs, bool flag_grad = false, bool flag_dgrad = false);
    void _calcObs(const double * vp, double &f, double &df, double * grad_E = nullptr, double * dgrad_E = nullptr);
public:
    StochasticReconfigurationTargetFunction(VMC &vmc, int64_t E_Nmc, int64_t grad_E_Nmc, bool useGradErr, double lambda_reg = 0.):
            nfm::NoisyFunctionWithGradient(vmc.getNVP(), useGradErr), _vmc(vmc),
            _E_Nmc(E_Nmc), _grad_E_Nmc(grad_E_Nmc), _lambda_reg(lambda_reg) {}

    ~StochasticReconfigurationTargetFunction() final = default;

    // NoisyFunctionWithGradient implementation
    nfm::NoisyValue f(const std::vector<double> &vp) final;
    void grad(const std::vector<double> &vp, nfm::NoisyGradient &grad) final;
    nfm::NoisyValue fgrad(const std::vector<double> &vp, nfm::NoisyGradient &grad) final;
};
} // namespace vmc

#endif
