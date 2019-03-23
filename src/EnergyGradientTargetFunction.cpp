#include "vmc/EnergyGradientTargetFunction.hpp"

#include "vmc/EnergyGradientMCObservable.hpp"
#include "vmc/MPIVMC.hpp"

void EnergyGradientTargetFunction::f(const double *vp, double &f, double &df)
{
    // set the variational parameters given as input
    _wf->setVP(vp);
    // perform the integral and store the values
    double obs[4];
    double dobs[4];
    MPIVMC::Integrate(_mci, _E_Nmc, obs, dobs, true, true);
    f = obs[0];
    df = dobs[0];

    if (_lambda_reg > 0.) { // compute the regularization term
        double norm = 0.;
        for (int i=0; i<_wf->getNVP(); ++i) { norm += vp[i]*vp[i]; }
        f += _lambda_reg*norm/_wf->getNVP();
    }
}


void EnergyGradientTargetFunction::grad(const double *vp, double *grad_E, double *dgrad_E)
{
    double f, df;
    fgrad(vp, f, df, grad_E, dgrad_E);
}


void EnergyGradientTargetFunction::fgrad(const double *vp, double &f, double &df, double *grad_E, double *dgrad_E)
{
    // set the variational parameters given as input
    _wf->setVP(vp);
    // add gradient obs to MCI
    _mci->addObservable( EnergyGradientMCObservable(_wf, _H), 1, 1, false, true ); // skipping equlibiration for gradients
    // perform the integral and store the values
    double obs[4 + 2*_wf->getNVP()];
    double dobs[4 + 2*_wf->getNVP()];
    MPIVMC::Integrate(_mci, _grad_E_Nmc, obs, dobs, true, true);
    // create pointers for ease of use and readability
    const double * const H = obs;
    const double * const dH = dobs;
    const double * const Oi = obs+4;
    const double * const dOi = dobs+4;
    const double * const HOi = obs+4+_wf->getNVP();
    const double * const dHOi = dobs+4+_wf->getNVP();
    f = H[0];
    df = dH[0];
    // compute direction (or gradient) to follow
    for (int i=0; i<_wf->getNVP(); ++i){
        grad_E[i] = 2.*( HOi[i] - H[0]*Oi[i] );
        dgrad_E[i] = 2.*( dHOi[i] + fabs(H[0]*Oi[i])*(dH[0]/H[0]+dOi[i]/Oi[i]) );
    }

    if (_lambda_reg > 0.) { // compute the regularization derivative
        const double fac = 2.*_lambda_reg/_wf->getNVP();
        for (int i=0; i<_wf->getNVP(); ++i) { grad_E[i] += fac*vp[i]; }
    }
    _mci->popObservable(); // remove the gradient obs (it will be deleted)
}
