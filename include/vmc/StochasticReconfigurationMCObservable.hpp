#ifndef VMC_STOCHASTICRECONFIGURATIONMCOBSERVABLE_HPP
#define VMC_STOCHASTICRECONFIGURATIONMCOBSERVABLE_HPP

#include "mci/ObservableFunctionInterface.hpp"
#include "mci/DependentObservableInterface.hpp"
#include "vmc/Hamiltonian.hpp"
#include "vmc/WaveFunction.hpp"
#include "vmc/DependencyHelpers.hpp"

namespace vmc
{
// MC Observable used to sample required ingredients for Stochastic Reconfiguration Gradients
//
// NOTE: Remember the necessary dependency binding just as in the case of Hamiltonian.
class StochasticReconfigurationMCObservable: public mci::ObservableFunctionInterface, public mci::DependentObservableInterface
{
protected:
    const int _nvp; // number of variational parameters (nobs will be twice of that)

    // These must be bound via registerDeps() (called by MCI) or bind methods
    const WaveFunction * _wf = nullptr; // WaveFunction to read derivatives from
    const double * _E = nullptr; // if set, should point to an array of len 4 where the energies are read from

    mci::ObservableFunctionInterface * _clone() const final
    {
        return new StochasticReconfigurationMCObservable(_ndim, _nvp);
    }
public:
    StochasticReconfigurationMCObservable(int ntotaldim, int nvp):
            mci::ObservableFunctionInterface(ntotaldim, 2*nvp + nvp*nvp, false),
            mci::DependentObservableInterface(true), _nvp(nvp) {}

    ~StochasticReconfigurationMCObservable() final = default;


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
        _E = fetchEnergyDep<Hamiltonian>(accus, selfIdx, "StochasticReconfigurationMCObservable::registerDeps");
        _wf = fetchWaveFunctionDep<WaveFunction>(pdfcont, "StochasticReconfigurationMCObservable::registerDeps");
    }

    void deregisterDeps() final
    {
        _E = nullptr;
        _wf = nullptr;
    }

    // mci::ObservableFunctionInterface implementation
    void observableFunction(const double * in, double * out) final
    {
        // out is made in this way (nvp is the number of variational parameters):
        // out[0:nvp-1] = Oi
        // out[nvp:2*nvp-1] = HOi
        // out[2*nvp:2*nvp+nvp*nvp-1] = OiOj    according to the rule OiOj[i, j] = OiOj[j + i*nvp]

        // local energy
        const double Hloc = _E[ElocID::ETot]; // use enum integer to get total energy

        // store the elements Oi and HOi
        for (int i = 0; i < _nvp; ++i) {
            out[i] = _wf->getVD1DivByWF(i); //  Oi
            out[i + _nvp] = Hloc*out[i];    // HOi
        }
        // store the elements OiOj
        for (int i = 0; i < _nvp; ++i) {
            for (int j = 0; j < _nvp; ++j) {
                out[2*_nvp + i*_nvp + j] = out[i]*out[j]; // OiOj
            }
        }
    }
};
} // namespace vmc

#endif
