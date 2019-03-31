#ifndef VMC_ENERGYGRADIENTTARGETFUNCTION_HPP
#define VMC_ENERGYGRADIENTTARGETFUNCTION_HPP

#include "vmc/Hamiltonian.hpp"
#include "vmc/WaveFunction.hpp"

#include "mci/MCIntegrator.hpp"
#include "nfm/NoisyFunction.hpp"

namespace vmc
{

class EnergyGradientTargetFunction: public nfm::NoisyFunctionWithGradient
{
protected:
    WaveFunction * const _wf;
    Hamiltonian * const _H;
    mci::MCI * const _mci;
    const int _E_Nmc;
    const int _grad_E_Nmc;
    const double _lambda_reg;

public:
    EnergyGradientTargetFunction(WaveFunction * wf, Hamiltonian * H, const int &E_Nmc, const int &grad_E_Nmc, mci::MCI * mci, const double &lambda_reg = 0.):
            nfm::NoisyFunctionWithGradient(wf->getNVP(), true),
            _wf(wf), _H(H), _mci(mci), _E_Nmc(E_Nmc), _grad_E_Nmc(grad_E_Nmc), _lambda_reg(lambda_reg) {}

    ~EnergyGradientTargetFunction() override = default;


    // NoisyFunctionWithGradient implementation
    nfm::NoisyValue f(const std::vector<double> &vp) override;
    void grad(const std::vector<double> &vp, std::vector<nfm::NoisyValue> &grad) override;
    nfm::NoisyValue fgrad(const std::vector<double> &vp, std::vector<nfm::NoisyValue> &grad) override;
};
} // namespace vmc

#endif
