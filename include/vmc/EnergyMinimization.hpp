#ifndef VMC_ENERGYMINIMIZATION_HPP
#define VMC_ENERGYMINIMIZATION_HPP

#include "vmc/VMC.hpp"
#include "vmc/EnergyGradientTargetFunction.hpp" // we use that as default
#include "nfm/NoisyFunMin.hpp"

namespace vmc
{

// Generic gradient-based minimizeEnergy function
//
// Pass vmc to be optimized and existing noisy optimizer NFM
// Set the Gradient Type to be used via template <> and pass
// the required constructor arguments.
template <typename GradT = EnergyGradientTargetFunction>
void minimizeEnergy(VMC &vmc, nfm::NFM &nfm, int E_NMC, int grad_E_NMC, bool useGradErr = true, double lambda_reg = 0.)
{
    GradT gradfun(vmc, E_NMC, grad_E_NMC, useGradErr, lambda_reg); // create gradient target function of type GradT with passed arguments
    std::vector<double> x0(static_cast<size_t>(vmc.getNVP()));
    vmc.getVP(x0.data()); // set x0 from wf VP
    nfm.findMin(gradfun, x0); // minimize energy
}
} // namespace vmc


#endif
