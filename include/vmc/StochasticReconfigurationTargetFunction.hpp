#ifndef VMC_STOCHASTICRECONFIGURATIONTARGETFUNCTION_HPP
#define VMC_STOCHASTICRECONFIGURATIONTARGETFUNCTION_HPP

#include "vmc/Hamiltonian.hpp"
#include "vmc/WaveFunction.hpp"

#include "mci/MCIntegrator.hpp"
#include "nfm/NoisyFunction.hpp"

#include <iostream>

namespace vmc
{

class StochasticReconfigurationTargetFunction: public NoisyFunctionWithGradient
{
protected:
    WaveFunction * const _wf;
    Hamiltonian * const _H;
    mci::MCI * const _mci;
    const int _Nmc;
    const double _lambda_reg;
    const bool _calcDGrad; // allows to disable calculation of gradient error

    void _integrate(const double * vp, double * obs, double * dobs, bool flag_grad = false, bool flag_dgrad = false);
    void _calcObs(const double * vp, double &f, double &df, double * grad_E = nullptr, double * dgrad_E = nullptr);
public:
    StochasticReconfigurationTargetFunction(WaveFunction * wf, Hamiltonian * H, mci::MCI * mci, const int &Nmc, const double &lambda_reg = 0., const bool calcDGrad = false):
        NoisyFunctionWithGradient(wf->getNVP()), _wf(wf), _H(H), _mci(mci), _Nmc(Nmc), _lambda_reg(lambda_reg), _calcDGrad(calcDGrad) {}

    ~StochasticReconfigurationTargetFunction() override= default;

    // NoisyFunctionWithGradient implementation
    void f(const double *vp, double &f, double &df) override;
    void grad(const double *vp, double *grad_E, double *dgrad_E) override;
    void fgrad(const double *vp, double &f, double &df, double *grad_E, double *dgrad_E) override;
};
} // namespace vmc

#endif
