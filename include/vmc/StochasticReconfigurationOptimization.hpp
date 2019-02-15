#ifndef STOCHASTIC_RECONFIGURATION_OPTIMIZATION
#define STOCHASTIC_RECONFIGURATION_OPTIMIZATION

#include "vmc/WFOptimization.hpp"
#include "mci/MCIntegrator.hpp"

#include <gsl/gsl_vector.h>

class StochasticReconfigurationOptimization: public WFOptimization{

private:
    const long _Nmc;
    const double _lambda_reg;
    const double _stepSize;

public:
    StochasticReconfigurationOptimization(WaveFunction * wf, Hamiltonian * H, const long &Nmc, MCI * mci, const double stepSize = 1., const double &lambda_reg = 0)
        : WFOptimization(wf, H, mci), _Nmc(Nmc), _lambda_reg(lambda_reg), _stepSize(stepSize) {}
    virtual ~StochasticReconfigurationOptimization(){}

    long getNmc(){return _Nmc;}
    double getLambdaReg(){return _lambda_reg;}

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
        double lambda_reg; // L2 regularization lambda

        void initFromOptimizer(StochasticReconfigurationOptimization * wfopt);
    };

    void fval(vmc_workspace &w, const double *vp, double &f, double &df);
    void grad(vmc_workspace &w, const double *vp, double *grad_E, double *dgrad_E = NULL);
    void fgrad(vmc_workspace &w, const double *vp, double &f, double &df, double *grad_E, double *dgrad_E = NULL);
};

#endif
