#ifndef STOCHASTIC_RECONFIGURATION_OPTIMIZATION
#define STOCHASTIC_RECONFIGURATION_OPTIMIZATION

#include "WFOptimization.hpp"
#include "MCIntegrator.hpp"

#include <gsl/gsl_vector.h>

class StochasticReconfigurationOptimization: public WFOptimization{

private:
    long _Nmc;

public:
    StochasticReconfigurationOptimization(WaveFunction * wf, Hamiltonian * H, const long &Nmc, MCI * mci): WFOptimization(wf, H, mci){
        _Nmc = Nmc;
    }
    virtual ~StochasticReconfigurationOptimization(){}

    long getNmc(){return _Nmc;}

    // optimization
    void optimizeWF();
};

namespace sropt_details {
    struct vmc_workspace
    {
        WaveFunction * wf;
        Hamiltonian * H;
        MCI * mci;
        long Nmc;

        void initFromOptimizer(StochasticReconfigurationOptimization * wfopt);
    };

    void fval(vmc_workspace &w, const double *vp, double &f, double &df);
    void grad(vmc_workspace &w, const double *vp, double *grad_E, double *dgrad_E = NULL);
    void fgrad(vmc_workspace &w, const double *vp, double &f, double &df, double *grad_E, double *dgrad_E = NULL);
};

#endif
