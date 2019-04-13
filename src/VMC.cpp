#include "vmc/VMC.hpp"
#include "vmc/MPIVMC.hpp"

namespace vmc
{

// --- derivative callback

VMC::DerivativeCallback::DerivativeCallback(WaveFunction * wf):
        mci::CallBackOnMoveInterface(wf->getNDim()), _wf(wf) {}

void VMC::DerivativeCallback::callBackFunction(const mci::WalkerState &wlk)
{
    if (wlk.accepted && wlk.needsObs) {
        _wf->computeAllDerivatives(wlk.xnew);
    }
}


// --- compute quantities

void VMC::computeEnergy(const int Nmc, double E[], double dE[], const bool doFindMRT2step, const bool doDecorrelation)
{
    MPIVMC::Integrate(_mci, Nmc, E, dE, doFindMRT2step, doDecorrelation);
}


// --- Optimization methods

void VMC::conjugateGradientOptimization(const int E_Nmc, const int grad_E_Nmc)
{
    ConjugateGradientOptimization opt(_wf, _H, E_Nmc, grad_E_Nmc, &_mci);
    opt.optimizeWF();
};

void VMC::stochasticReconfigurationOptimization(const int Nmc, const double stepSize, const bool flag_dgrad)
{
    StochasticReconfigurationOptimization opt(_wf, _H, Nmc, &_mci, stepSize, flag_dgrad);
    opt.optimizeWF();
};

void VMC::adamOptimization(const int Nmc, const bool useSR, const bool useGradientError, const int max_n_const_values, const bool useAveraging,
                           const double lambda, const double alpha, const double beta1, const double beta2, const double epsilon)
{
    AdamOptimization opt(_wf, _H, &_mci, Nmc, useSR, useGradientError, max_n_const_values, useAveraging, lambda, alpha, beta1, beta2, epsilon);
    opt.optimizeWF();
};


void VMC::simulatedAnnealingOptimization(const int Nmc, const double iota, const double kappa, const double lambda, gsl_siman_params_t &params)
{
    SimulatedAnnealingOptimization opt(_wf, _H, Nmc, &_mci, iota, kappa, lambda, params);
    opt.optimizeWF();
};

void VMC::nmsimplexOptimization(const int Nmc, const double iota, const double kappa, const double lambda, const double rstart, const double rend, const size_t max_n_iter)
{
    NMSimplexOptimization opt(_wf, _H, &_mci, Nmc, iota, kappa, lambda, rstart, rend, max_n_iter);
    opt.optimizeWF();
};
} // namespace vmc