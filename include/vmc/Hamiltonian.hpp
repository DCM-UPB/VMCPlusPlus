#ifndef VMC_HAMILTONIAN_HPP
#define VMC_HAMILTONIAN_HPP

#include "mci/ObservableFunctionInterface.hpp"
#include "mci/DependentObservableInterface.hpp"
#include "vmc/WaveFunction.hpp"
#include "vmc/DependencyHelpers.hpp"

namespace vmc
{

// Base class for Hamiltonian MCI observables
//
// To finalize the Hamiltonian you only need to implement the localPotentialEnergy() method (and _clone()).
//
// NOTE:
// Hamiltonians are strongly coupled with the wave function / sampling function and therefore only
// work when a WaveFunction is known. However, requiring a permanent binding on construction is
// undesirable in our context. Instead we rely on registration mechanisms of MCI.
// Remember that this means calling the observableFunction() and kinetic energy methods is not safe
// unless Hamiltonian is known to be bound or checked with isBound(). For use of the Hamiltonian outside
// of MCI you may use the bindWaveFunction method to bind a WaveFunction.
//

// enum to use if you don't want to remember the energy index mapping of MCI observable array (see observableFunction())
enum ElocID { ETot = 0, EPot = 1, EKinPB = 2, EKinJF = 3 };

// Hamiltonian: obs[0]=Totale Energy, obs[1]=Potential Energy, obs[2]=Kinetic Energy (PB), obs[3]=Kinetic Energy (JF)
class Hamiltonian: public mci::ObservableFunctionInterface, public mci::DependentObservableInterface
{
protected:
    const int _nspacedim;
    const int _npart;

    const bool _flag_PBKE;

    const WaveFunction * _wf = nullptr; // pointer to the wf we depend on

    Hamiltonian(int nspacedim, int npart, bool usePBKE = true /* use JF+PB KE or only JF */):
            mci::ObservableFunctionInterface(nspacedim*npart, 4, false), mci::DependentObservableInterface(true),
            _nspacedim(nspacedim), _npart(npart), _flag_PBKE(usePBKE) {}

public:
    ~Hamiltonian() override = default;

    int getNSpaceDim() const { return _nspacedim; }
    int getTotalNDim() const { return getNDim(); }
    int getNPart() const { return _npart; }

    bool hasPBKE() const { return _flag_PBKE; }

    // use this if you need to check whether observable is bound (i.e. can be fully used)
    bool isBound() const { return (_wf != nullptr); }

    // Simple method to bind a wave function (e.g. for testing)
    void bindWaveFunction(WaveFunction * wf) { _wf = wf; } // if wf==nullptr, isBound() will be false afterwards

    // Methods to register/deregister WaveFunction dependency (called by MCI)
    void registerDeps(const mci::SamplingFunctionContainer &pdfcont, const std::vector<mci::AccumulatorInterface *> &/*accus*/, int/*selfIdx*/) override
    {
        // use helper
        _wf = fetchWaveFunctionDep<WaveFunction>(pdfcont, "Hamiltonian::registerDeps");
    }

    void deregisterDeps() override
    {
        _wf = nullptr;
    }


    // --- Energy methods

    // Potential energy --- MUST BE IMPLEMENTED
    virtual double localPotentialEnergy(const double * r) = 0;

    double localPBKineticEnergy(const double * /*r*/)
    {
        double ekin = 0.;
        for (int i = 0; i < _ndim; ++i) {
            ekin += _wf->getD2DivByWF(i);
        }
        return (-0.5*ekin);
    }

    double localJFKineticEnergy(const double * /*r*/)
    {
        double ekin = 0.;
        for (int i = 0; i < _ndim; ++i) {
            const double foo = _wf->getD1DivByWF(i);
            ekin += foo*foo;
        }
        return (0.5*ekin);
    }

    void observableFunction(const double * in, double * out) override
    {
        out[ElocID::EKinJF] = localJFKineticEnergy(in);
        out[ElocID::EKinPB] = _flag_PBKE ? localPBKineticEnergy(in) : 0.;
        out[ElocID::EPot] = localPotentialEnergy(in);
        out[ElocID::ETot] = _flag_PBKE ? out[ElocID::EPot] + out[EKinPB] : out[ElocID::EPot] + out[ElocID::EKinJF];
    }
};
} // namespace vmc

#endif
