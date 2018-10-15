#ifndef NOISY_STOCHASTIC_RECONFIGURATION_TARGET_FUNCTION
#define NOISY_STOCHASTIC_RECONFIGURATION_TARGET_FUNCTION

#include "StochasticReconfigurationOptimization.hpp"
#include "MCIntegrator.hpp"
#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "NoisyFunction.hpp"

#include <iostream>



class NoisyStochasticReconfigurationTargetFunction: public NoisyFunctionWithGradient
{
protected:
    WaveFunction * _wf;
    Hamiltonian * _H;
    MCI * _mci;
    const long _Nmc;
    const double _lambda_reg;
    const bool _calcDGrad; // allows to disable calculation of gradient error

public:
    NoisyStochasticReconfigurationTargetFunction(WaveFunction * wf, Hamiltonian * H, MCI * mci, const long &Nmc, const double &lambda_reg = 0., const bool calcDGrad = false):
        NoisyFunctionWithGradient(wf->getNVP()), _Nmc(Nmc), _lambda_reg(lambda_reg), _calcDGrad(calcDGrad) {
        _wf = wf;
        _H = H;
        _mci = mci;
    }

    virtual ~NoisyStochasticReconfigurationTargetFunction(){}

    // NoisyFunctionWithGradient implementation
    void f(const double *vp, double &f, double &df)
    {
        sropt_details::vmc_workspace w;
        w.wf = _wf;
        w.H = _H;
        w.mci = _mci;
        w.Nmc = _Nmc;
        w.lambda_reg = _lambda_reg;
        sropt_details::fval(w, vp, f, df);

    }

    void grad(const double *vp, double *grad_E, double *dgrad_E)
    {
        sropt_details::vmc_workspace w;
        w.wf = _wf;
        w.H = _H;
        w.mci = _mci;
        w.Nmc = _Nmc;
        w.lambda_reg = _lambda_reg;
        sropt_details::grad(w, vp, grad_E, _calcDGrad ? dgrad_E : NULL);
    }

    void fgrad(const double *vp, double &f, double &df, double *grad_E, double *dgrad_E)
    {
        sropt_details::vmc_workspace w;
        w.wf = _wf;
        w.H = _H;
        w.mci = _mci;
        w.Nmc = _Nmc;
        w.lambda_reg = _lambda_reg;
        sropt_details::fgrad(w, vp, f, df, grad_E, _calcDGrad ? dgrad_E : NULL);
    }

};


#endif
