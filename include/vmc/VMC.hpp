#ifndef VMC_VMC_HPP
#define VMC_VMC_HPP

#include "vmc/Hamiltonian.hpp"
#include "vmc/WaveFunction.hpp"
#include "mci/MCIntegrator.hpp"

#include <stdexcept>
#include <memory>

namespace vmc
{

class VMC
{
protected:
    mci::MCI _mci; // the contained MC integrator

    // Internal pointers to wf and H inside MCI
    WaveFunction * const _wf;
    Hamiltonian * const _H;

public:
    // Constructors
    VMC(std::unique_ptr<WaveFunction> wf, std::unique_ptr<Hamiltonian> H);  // move unique pointers into VMC
    VMC(const WaveFunction &wf, const Hamiltonian &H); // clone provided wf and H

    // You may directly access and edit the MCI object
    // NOTE: This allows for some unsafe operations. In particular. adding extra
    // wave functions as sampling functions or clearing/popping sampling functions,
    // observables or callbacks added by VMC is strictly prohibited.
    mci::MCI &getMCI() { return _mci; }

    // Late access to contained WF/H
    WaveFunction &getWF() const { return *_wf; }
    Hamiltonian &getH() const { return *_H; }

    // Other accessors
    int getNSpaceDim() const { return _H->getNSpaceDim(); }
    int getNParticles() const { return _H->getNPart(); }
    int getNTotalDim() const { return _H->getTotalNDim(); }
    int getNVP() const { return _wf->getNVP(); }

    void setVP(const double * vp) { _wf->setVP(vp); }
    void getVP(double * vp) const { _wf->getVP(vp); }


    // Computation of the energy according to contained Hamiltonian and WaveFunction
    // Other contained observables will be calculated as well and stored behind the energy values
    void computeEnergy(int Nmc, double * E, double * dE, bool doFindMRT2step = true, bool doDecorrelation = true);
};
} // namespace vmc

#endif
