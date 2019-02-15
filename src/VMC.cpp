#include "vmc/VMC.hpp"
#include "vmc/MPIVMC.hpp"

// --- compute quantities

void VMC::computeVariationalEnergy(const long & Nmc, double * E, double * dE, const bool doFindMRT2step, const bool doDecorrelation)
{
    getMCI()->clearSamplingFunctions(); getMCI()->addSamplingFunction(_wf);
    getMCI()->clearObservables(); getMCI()->addObservable(_H);
    MPIVMC::Integrate(getMCI(), Nmc, E, dE, doFindMRT2step, doDecorrelation);
}


// --- Optimization methods

void VMC::conjugateGradientOptimization(const long &E_Nmc, const long &grad_E_Nmc)
{
    ConjugateGradientOptimization * opt = new ConjugateGradientOptimization(_wf, _H, E_Nmc, grad_E_Nmc, getMCI());
    opt->optimizeWF();
    delete opt;
};

void VMC::stochasticReconfigurationOptimization(const long &Nmc, const double stepSize, const bool flag_noisy)
{
    WFOptimization * opt;
    if (flag_noisy) opt = new NoisyStochasticReconfigurationOptimization(_wf, _H, Nmc, getMCI(), stepSize);
    else opt = new StochasticReconfigurationOptimization(_wf, _H, Nmc, getMCI(), stepSize);
    opt->optimizeWF();
    delete opt;
};

void VMC::adamOptimization(const long &Nmc, const bool useSR, const bool useGradientError, const size_t &max_n_const_values, const bool useAveraging,
    const double &lambda, const double &alpha, const double &beta1, const double &beta2, const double &epsilon)
{
    AdamOptimization * opt = new AdamOptimization(_wf, _H, getMCI(), Nmc, useSR, useGradientError, max_n_const_values, useAveraging, lambda, alpha, beta1, beta2, epsilon);
    opt->optimizeWF();
    delete opt;
};


void VMC::simulatedAnnealingOptimization(const long &Nmc, const double &iota, const double &kappa, const double &lambda, gsl_siman_params_t &params)
{
    SimulatedAnnealingOptimization * opt = new SimulatedAnnealingOptimization(_wf, _H, Nmc, getMCI(), iota, kappa, lambda, params);
    opt->optimizeWF();
    delete opt;
};

void VMC::nmsimplexOptimization(const long &Nmc, const double &iota, const double &kappa, const double &lambda, const double &rstart, const double &rend, const size_t &max_n_iter)
{
    NMSimplexOptimization * opt = new NMSimplexOptimization(_wf, _H, _mci, Nmc, iota, kappa, lambda, rstart, rend, max_n_iter);
    opt->optimizeWF();
    delete opt;
};
