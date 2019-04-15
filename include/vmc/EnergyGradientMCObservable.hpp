#ifndef VMC_ENERGYGRADIENTMCOBSERVABLE_HPP
#define VMC_ENERGYGRADIENTMCOBSERVABLE_HPP

#include "mci/ObservableFunctionInterface.hpp"
#include "mci/DependentObservableInterface.hpp"
#include "vmc/Hamiltonian.hpp"
#include "vmc/WaveFunction.hpp"
#include "vmc/DependencyHelpers.hpp"

namespace vmc
{

// MC Observable used to sample required ingredients for standard Hellman-Feynman energy gradients
//
// NOTE: Remember the necessary dependency binding just as in the case of Hamiltonian.
class EnergyGradientMCObservable: public mci::ObservableFunctionInterface, public mci::DependentObservableInterface
{
protected:
    const int _nvp; // number of variational parameters (nobs will be twice of that)

    // These must be bound via registerDeps() (called by MCI) or bind methods
    const WaveFunction * _wf = nullptr; // WaveFunction to read derivatives from
    const double * _E = nullptr; // if set, should point to an array of len 4 where the energies are read from

    mci::ObservableFunctionInterface * _clone() const final
    {
        return new EnergyGradientMCObservable(_ndim, _nvp);
    }
public:
    EnergyGradientMCObservable(int ntotaldim, int nvp):
            mci::ObservableFunctionInterface(ntotaldim, 2*nvp, false),
            mci::DependentObservableInterface(true), _nvp(nvp) {}

    ~EnergyGradientMCObservable() final = default;

    // use this if you need to check whether required objects are bound (i.e. may be fully used)
    bool isBound() const { return (_E != nullptr && _wf != nullptr); }

    // Simple method to bind a wave function (e.g. for testing)
    void bindDependencies(const double E[], WaveFunction * wf) // if any is nullptr, isBound() will be false afterwards
    {
        _E = E;
        _wf = wf;
    }

    // Methods to register/deregister WaveFunction/Hamiltonian dependency (called by MCI)
    void registerDeps(const mci::SamplingFunctionContainer &pdfcont, const std::vector<mci::AccumulatorInterface *> &accus, int selfIdx) final
    {
        // use dependency helpers
        _E = fetchEnergyDep<Hamiltonian>(accus, selfIdx, "EnergyGradientMCObservable::registerDeps");
        _wf = fetchWaveFunctionDep<WaveFunction>(pdfcont, "EnergyGradientMCObservable::registerDeps");
    }

    void deregisterDeps() final
    {
        _E = nullptr;
        _wf = nullptr;
    }

    // mci::ObservableFunctionInterface implementation
    void observableFunction(const double * in, double * out) final
    {
        // obs[0 : wf->getNVP()-1] = Variational Derivative of the Wave Function
        // obs[wf->getNVP() : 2*wf->getNVP()-1] = Local Energy times the Variational Derivative of the Wave Function
        const double Hloc = _E[ElocID::ETot]; // use enum integer to get total energy
        for (int i = 0; i < _wf->getNVP(); ++i) {
            out[i] = _wf->getVD1DivByWF(i);
            out[i + _wf->getNVP()] = out[i]*Hloc;
        }
    }
};
} // namespace vmc

#endif
