#ifndef VMC_DEPENDENCYHELPERS_HPP
#define VMC_DEPENDENCYHELPERS_HPP

#include "mci/DependentObservableInterface.hpp"

#include <stdexcept>

namespace vmc
{

// Helpers functions for VMC's observables to safely fetch required WaveFunction or Hamiltonian dependencies
//
// The helpers are templated so that they can also be used for possible future extended interfaces.

// Safely get a valid pointer to the VMC WaveFunction in pdfcont
template <class WFType>
const WFType * fetchWaveFunctionDep(const mci::SamplingFunctionContainer &pdfcont, const std::string &callerID)
{
    if (pdfcont.size() != 1) {
        throw std::runtime_error("[" + callerID + "] The passed PDF container does not contain exactly one sampling function.");
    }
    auto wf = dynamic_cast<WFType *>(&pdfcont.getSamplingFunction(0));
    if (wf == nullptr) {
        throw std::runtime_error("[" + callerID + "] The passed sampling function was not derived from WaveFunction type.");
    }
    return wf;
}

// Safely get a valid pointer to the VMC Hamiltonian energy array (and check that the Hamiltonian is indeed of specified type)
template <class HType>
const double * fetchEnergyDep(const std::vector<mci::AccumulatorInterface *> &accus, int selfIdx, const std::string &callerID)
{
    if (accus.size() < 2) {
        throw std::runtime_error("[" + callerID + "] The passed obs container does not contain at least two observables (including this).");
    }
    // check that first observable is derived from Hamiltonian
    const auto H = dynamic_cast<HType *>(&accus[0]->getObservableFunction()); // just stored for a check
    if (H == nullptr) {
        throw std::runtime_error("[" + callerID + "] The first passed observable is not derived from Hamiltonian type.");
    }
    if (!mci::DependentObservableInterface::isObsDepValid(accus, selfIdx, 0)) {
        throw std::runtime_error("[" + callerID + "] The accumulators of Hamiltonian and this have unsynced nskip settings.");
    }
    return accus[0]->getObsValues(); // now we can be sure that everything is fine
}
} // namespace vmc

#endif
