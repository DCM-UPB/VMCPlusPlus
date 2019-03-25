#ifndef VMC_SIMULATEDANNEALINGOPTIMIZATION_HPP
#define VMC_SIMULATEDANNEALINGOPTIMIZATION_HPP

#include "vmc/WFOptimization.hpp"

#include <gsl/gsl_siman.h>

namespace vmc
{

class SimulatedAnnealingOptimization: public WFOptimization
{
public:
    SimulatedAnnealingOptimization(WaveFunction * wf, Hamiltonian * H, const int &Nmc, mci::MCI * mci, const double &iota, const double &kappa, const double &lambda, gsl_siman_params_t &params);

    ~SimulatedAnnealingOptimization() override= default;

    // optimization
    void optimizeWF() override;
};
} // namespace vmc

#endif
