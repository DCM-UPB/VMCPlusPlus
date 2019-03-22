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


class VMC
{
protected:
    WaveFunction * const _wf;
    Hamiltonian * const _H;
    mci::MCI * _mci;


public:
    VMC(WaveFunction * wf, Hamiltonian * H):
        _wf(wf), _H(H)
    {
        if (_wf->getTotalNDim() != _H->getTotalNDim()) {
            throw std::invalid_argument( "Error VMC: ndim different between wf and H" );
        }
        _mci = new mci::MCI(_H->getTotalNDim());
    }

    ~VMC(){
        delete _mci;
    }


    // Monte Carlo Integral within VMC should be performed using the MCI object provided by VMC
    mci::MCI * getMCI(){return _mci;}

    // Computation of the variational energy
    void computeVariationalEnergy(const int & Nmc, double * E, double * dE, bool doFindMRT2step = true, bool doDecorrelation = true);


    // Wave Function Optimization Methods
    void conjugateGradientOptimization(const int &E_Nmc, const int &grad_E_Nmc);

    void stochasticReconfigurationOptimization(const int &Nmc, double stepSize = 1., bool flag_dgrad = false); // calc&use gradient error?

    void adamOptimization(const int &Nmc, bool useSR = false, bool useGradientError = false, const size_t &max_n_const_values = 20, bool useAveraging = false,
        const double &lambda = 0, const double &alpha = 0.001, const double &beta1 = 0.9, const double &beta2 = 0.999, const double &epsilon = 10e-8);

    void simulatedAnnealingOptimization(const int &Nmc, const double &iota, const double &kappa, const double &lambda, gsl_siman_params_t &params);

    void nmsimplexOptimization(const int &Nmc, const double &iota, const double &kappa, const double &lambda, const double &rstart = 1.0, const double &rend = 0.01, const size_t &max_n_iter = 0);
};


#endif
