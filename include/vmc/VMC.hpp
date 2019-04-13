#ifndef VMC_VMC_HPP
#define VMC_VMC_HPP

#include "vmc/AdamOptimization.hpp"
#include "vmc/ConjugateGradientOptimization.hpp"
#include "vmc/Hamiltonian.hpp"
#include "vmc/NMSimplexOptimization.hpp"
#include "vmc/SimulatedAnnealingOptimization.hpp"
#include "vmc/StochasticReconfigurationOptimization.hpp"
#include "vmc/StochasticReconfigurationOptimization.hpp"
#include "vmc/WaveFunction.hpp"

#include "mci/MCIntegrator.hpp"

#include <gsl/gsl_siman.h>
#include <stdexcept>
#include <memory>

namespace vmc
{

class VMC
{
    class DerivativeCallback: public mci::CallBackOnMoveInterface
        // Small internal helper (MCI CallBack Function)
        // Triggers WF derivative computation after MC move is accepted
    {
    protected:
        WaveFunction * const _wf;

        mci::CallBackOnMoveInterface * _clone() const final
        {
            return new DerivativeCallback(_wf);
        }

    public:
        explicit DerivativeCallback(WaveFunction * wf);

        void callBackFunction(const mci::WalkerState &wlk) final;
    };

protected:
    WaveFunction * const _wf;
    Hamiltonian * const _H;
    mci::MCI _mci;

public:
    VMC(std::unique_ptr<WaveFunction> wf, std::unique_ptr<Hamiltonian> H): // move unique pointers into VMC
            _wf(wf.get()/*remains valid (until destruct)*/), _H(H.get()), _mci(H->getTotalNDim())
    {
        if (_wf->getTotalNDim() != _H->getTotalNDim()) {
            throw std::invalid_argument("[VMC] ndim different between wf and H");
        }
        _mci.addSamplingFunction(std::move(wf));
        _mci.addObservable(std::move(H));
        //_mci.addSamplingFunction(std::unique_ptr<mci::SamplingFunctionInterface>(wf));
        //_mci.addObservable(std::unique_ptr<mci::ObservableFunctionInterface>(H));
        _mci.addCallBack(std::make_unique<DerivativeCallback>(_wf));
    }

    VMC(const WaveFunction &wf, const Hamiltonian &H):
            _wf(dynamic_cast<WaveFunction *>(wf.clone().release())), /*clone returns SamplingFunction ptr*/
            _H(dynamic_cast<Hamiltonian *>(H.clone().release())), /*clone returns ObservableFunction ptr*/
            _mci(H.getTotalNDim())
    {
        if (_wf == nullptr) {
            throw std::runtime_error("[VMC] WaveFunction's clone() did not produce a type derived from WaveFunction.");
        }
        if (_H == nullptr) {
            throw std::runtime_error("[VMC] Hamiltonian's clone() did not produce a type derived from Hamiltonian.");
        }
        if (_wf->getTotalNDim() != _H->getTotalNDim()) {
            throw std::invalid_argument("Error VMC: ndim different between wf and H");
        }
        _mci.addSamplingFunction(std::unique_ptr<mci::SamplingFunctionInterface>(_wf));
        _mci.addObservable(std::unique_ptr<mci::ObservableFunctionInterface>(_H));
        _mci.addCallBack(std::make_unique<DerivativeCallback>(_wf));
    }

    /*~VMC()
    {
        while (_mci.getNPDF() > 1) {
            // in case more pfs than _wf were added (e.g. via getMCI())
            _mci.popSamplingFunction();
        }
        auto wf = _mci.popSamplingFunction(); // reacquire unique pointer

        while (_mci.getNObs() > 1) {
            // in case something more than _wf was added (e.g. via getMCI())
            _mci.popObservable();
        }
        auto obs = _mci.popObservable(); // reacquire unique pointer
    }*/


    // You may directly access and edit the MCI object
    // NOTE: This allows for some unsafe operations. In particular. adding extra
    // wave functions as sampling functions or clearing/popping sampling functions,
    // observables or callbacks added by VMC is strictly prohibited.
    mci::MCI &getMCI() { return _mci; }

    // Late access to contained WF/H
    WaveFunction &getWF() const { return *_wf; }
    Hamiltonian &getH() const { return *_H; }

    // Computation of the energy
    void computeEnergy(int Nmc, double E[], double dE[], bool doFindMRT2step = true, bool doDecorrelation = true);


    // Wave Function Optimization Methods
    void conjugateGradientOptimization(int E_Nmc, int grad_E_Nmc);

    void stochasticReconfigurationOptimization(int Nmc, double stepSize = 1., bool flag_dgrad = false); // calc&use gradient error?

    void adamOptimization(int Nmc, bool useSR = false, bool useGradientError = false, int max_n_const_values = 20, bool useAveraging = false,
                          double lambda = 0, double alpha = 0.001, double beta1 = 0.9, double beta2 = 0.999, double epsilon = 10e-8);

    void simulatedAnnealingOptimization(int Nmc, double iota, double kappa, double lambda, gsl_siman_params_t &params);

    void nmsimplexOptimization(int Nmc, double iota, double kappa, double lambda, double rstart = 1.0, double rend = 0.01, size_t max_n_iter = 0);
};
} // namespace vmc

#endif
