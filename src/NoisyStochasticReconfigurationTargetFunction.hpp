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
    long _Nmc;
    MCI * _mci;

public:
    NoisyStochasticReconfigurationTargetFunction(WaveFunction * wf, Hamiltonian * H, const long & Nmc, MCI * mci):
        NoisyFunctionWithGradient(wf->getNVP()){
        _wf = wf;
        _H = H;
        _Nmc = Nmc;
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
        sropt_details::fval(w, vp, f, df);

    }

    void grad(const double *vp, double *grad_E, double *dgrad_E)
    {
        sropt_details::vmc_workspace w;
        w.wf = _wf;
        w.H = _H;
        w.mci = _mci;
        w.Nmc = _Nmc;
        sropt_details::grad(w, vp, grad_E, dgrad_E);
    }
};


#endif
