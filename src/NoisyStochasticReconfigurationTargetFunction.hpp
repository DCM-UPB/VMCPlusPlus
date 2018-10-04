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
    long _Nmc, _grad_E_Nmc;
    MCI * _mci;
    bool _calcDGrad; // allows to disable calculation of gradient error

public:
    NoisyStochasticReconfigurationTargetFunction(WaveFunction * wf, Hamiltonian * H, MCI * mci, const long &Nmc, const long &grad_E_Nmc = 0, const bool calcDGrad = true):
        NoisyFunctionWithGradient(wf->getNVP()){
        _wf = wf;
        _H = H;
        _mci = mci;
        _Nmc = Nmc;
        _grad_E_Nmc = (grad_E_Nmc > 0) ? grad_E_Nmc : Nmc; // use same number of steps for gradient if no extra number provided
        _calcDGrad = calcDGrad;
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
        sropt_details::fval(w, vp, f, df);

    }

    void grad(const double *vp, double *grad_E, double *dgrad_E)
    {
        sropt_details::vmc_workspace w;
        w.wf = _wf;
        w.H = _H;
        w.mci = _mci;
        w.Nmc = _grad_E_Nmc;
        sropt_details::grad(w, vp, grad_E, _calcDGrad ? dgrad_E : NULL);
    }
};


#endif
