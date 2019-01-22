#ifndef CONJUGATE_GRADIENT_TARGET_FUNCTION
#define CONJUGATE_GRADIENT_TARGET_FUNCTION


#include "ConjugateGradientMCObservable.hpp"
#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"

#include "MCIntegrator.hpp"
#include "NoisyFunction.hpp"
#include "MPIVMC.hpp"

class ConjugateGradientTargetFunction: public NoisyFunctionWithGradient
{
protected:
    WaveFunction * _wf;
    Hamiltonian * _H;
    MCI * _mci;
    const long _E_Nmc;
    const long _grad_E_Nmc;
    const double _lambda_reg;

public:
    ConjugateGradientTargetFunction(WaveFunction * wf, Hamiltonian * H, const long & E_Nmc, const long &grad_E_Nmc, MCI * mci, const double &lambda_reg = 0.):
        NoisyFunctionWithGradient(wf->getNVP()), _E_Nmc(E_Nmc), _grad_E_Nmc(grad_E_Nmc), _lambda_reg(lambda_reg) {
        _wf = wf;
        _H = H;
        _mci = mci;
    }

    virtual ~ConjugateGradientTargetFunction(){}


    // NoisyFunctionWithGradient implementation
    void f(const double *vp, double &f, double &df){
        // set the variational parameters given as input
        _wf->setVP(vp);
        // set up the MC integrator
        _mci->clearSamplingFunctions(); _mci->addSamplingFunction(_wf);
        _mci->clearObservables(); _mci->addObservable(_H);
        // perform the integral and store the values
        double * obs = new double[4];
        double * dobs = new double[4];
        MPIVMC::Integrate(_mci, _E_Nmc, obs, dobs, true, true);
        f = obs[0];
        df = dobs[0];

        if (_lambda_reg > 0.) { // compute the regularization term
            double norm = 0.;
            for (int i=0; i<_wf->getNVP(); ++i) norm += vp[i]*vp[i];
            f += _lambda_reg*norm/_wf->getNVP();
        }

        // free resources
        delete [] dobs;
        delete [] obs;
    }

    void grad(const double *vp, double *grad_E, double *dgrad_E)
    {
        double f, df;
        fgrad(vp, f, df, grad_E, dgrad_E);
    }

    void fgrad(const double *vp, double &f, double &df, double *grad_E, double *dgrad_E){
        // set the variational parameters given as input
        _wf->setVP(vp);
        // set up the MC integrator
        _mci->clearSamplingFunctions(); _mci->addSamplingFunction(_wf);
        _mci->clearObservables(); _mci->addObservable(_H);
        ConjugateGradientMCObservable * mc_obs = new ConjugateGradientMCObservable(_wf, _H);
        _mci->addObservable(mc_obs);
        // perform the integral and store the values
        double * obs = new double[4 + 2*_wf->getNVP()];
        double * dobs = new double[4 + 2*_wf->getNVP()];
        MPIVMC::Integrate(_mci, _grad_E_Nmc, obs, dobs, true, false); // skipping decorrelation for gradients
        // create pointers for ease of use and readability
        double * H = obs;
        double * dH = dobs;
        double * Oi = obs+4;
        double * dOi = dobs+4;
        double * HOi = obs+4+_wf->getNVP();
        double * dHOi = dobs+4+_wf->getNVP();
        f = H[0];
        df = dH[0];
        // compute direction (or gradient) to follow
        for (int i=0; i<_wf->getNVP(); ++i){
            grad_E[i] = 2.*( HOi[i] - H[0]*Oi[i] );
            dgrad_E[i] = 2.*( dHOi[i] + fabs(H[0]*Oi[i])*(dH[0]/H[0]+dOi[i]/Oi[i]) );
        }

        if (_lambda_reg > 0.) { // compute the regularization derivative
            const double fac = 2.*_lambda_reg/_wf->getNVP();
            for (int i=0; i<_wf->getNVP(); ++i) grad_E[i] += fac*vp[i];
        }

        // free resources
        delete mc_obs;
        delete[] dobs;
        delete[] obs;
    }
};


#endif
