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


    // Wave Function Optimization Methods
    /*  void conjugateGradientOptimization(int E_Nmc, int grad_E_Nmc);

      void stochasticReconfigurationOptimization(int Nmc, double stepSize = 1., bool flag_dgrad = false); // calc&use gradient error?

      void adamOptimization(int Nmc, bool useSR = false, bool useGradientError = false, int max_n_const_values = 20, bool useAveraging = false,
                            double lambda = 0, double alpha = 0.001, double beta1 = 0.9, double beta2 = 0.999, double epsilon = 10e-8);

      void simulatedAnnealingOptimization(int Nmc, double iota, double kappa, double lambda, gsl_siman_params_t &params);

      void nmsimplexOptimization(int Nmc, double iota, double kappa, double lambda, double rstart = 1.0, double rend = 0.01, size_t max_n_iter = 0);*/
};
} // namespace vmc

#endif
