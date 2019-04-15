#ifndef VMC_SIMULATEDANNEALINGMINIMIZATION_HPP
#define VMC_SIMULATEDANNEALINGMINIMIZATION_HPP

#include "vmc/VMC.hpp"

#include <gsl/gsl_siman.h>

namespace vmc
{

class SimulatedAnnealingMinimization
{
public:
    SimulatedAnnealingMinimization(int Nmc, double iota, double kappa, double lambda, gsl_siman_params_t &params);

    // optimization
    void minimizeEnergy(VMC &vmc);
};
} // namespace vmc

#endif
