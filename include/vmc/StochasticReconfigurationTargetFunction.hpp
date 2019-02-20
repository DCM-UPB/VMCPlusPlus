#ifndef STOCHASTIC_RECONFIGURATION_TARGET_FUNCTION
#define STOCHASTIC_RECONFIGURATION_TARGET_FUNCTION

#include "vmc/WaveFunction.hpp"
#include "vmc/Hamiltonian.hpp"

#include "nfm/NoisyFunction.hpp"
#include "mci/MCIntegrator.hpp"

#include <iostream>



class StochasticReconfigurationTargetFunction: public NoisyFunctionWithGradient
{
protected:
    WaveFunction * const _wf;
    Hamiltonian * const _H;
    MCI * const _mci;
    const long _Nmc;
    const double _lambda_reg;
    const bool _calcDGrad; // allows to disable calculation of gradient error

    void _integrate(const double * const vp, double * const obs, double * const dobs, const bool flag_grad = false);
    void _calcObs(const double * const vp, double &f, double &df, double * const grad_E = NULL, double * const dgrad_E = NULL);
public:
    StochasticReconfigurationTargetFunction(WaveFunction * wf, Hamiltonian * H, MCI * mci, const long &Nmc, const double &lambda_reg = 0., const bool calcDGrad = false):
        NoisyFunctionWithGradient(wf->getNVP()), _wf(wf), _H(H), _mci(mci), _Nmc(Nmc), _lambda_reg(lambda_reg), _calcDGrad(calcDGrad) {}

    virtual ~StochasticReconfigurationTargetFunction(){}

    // NoisyFunctionWithGradient implementation
    void f(const double *vp, double &f, double &df);
    void grad(const double *vp, double *grad_E, double *dgrad_E);
    void fgrad(const double *vp, double &f, double &df, double *grad_E, double *dgrad_E);
};


#endif
