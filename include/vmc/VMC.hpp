#ifndef VMC_CLASS
#define VMC_CLASS

#include "vmc/WaveFunction.hpp"
#include "vmc/Hamiltonian.hpp"
#include "vmc/ConjugateGradientOptimization.hpp"
#include "vmc/StochasticReconfigurationOptimization.hpp"
#include "vmc/StochasticReconfigurationOptimization.hpp"
#include "vmc/AdamOptimization.hpp"
#include "vmc/SimulatedAnnealingOptimization.hpp"
#include "vmc/NMSimplexOptimization.hpp"

#include "mci/MCIntegrator.hpp"

#include <stdexcept>
#include <gsl/gsl_siman.h>


class VMC{
protected:
    WaveFunction * _wf;
    Hamiltonian * _H;
    MCI * _mci;


public:
    VMC(WaveFunction * wf, Hamiltonian * H){
        using namespace std;
        if (wf->getTotalNDim() != H->getTotalNDim())
            throw std::invalid_argument( "Error VMC: ndim different between wf and H" );
        _wf=wf;
        _H=H;
        _mci = new MCI(_H->getTotalNDim());
        _mci->addCallBackOnAcceptance(_wf);
    }

    ~VMC(){
        delete _mci;
    }


    // Monte Carlo Integral within VMC should be performed using the MCI object provided by VMC
    MCI * getMCI(){return _mci;}


    // Computation of the variational energy
    void computeVariationalEnergy(const long & Nmc, double * E, double * dE, const bool doFindMRT2step = true, const bool doDecorrelation = true);


    // Wave Function Optimization Methods
    void conjugateGradientOptimization(const long &E_Nmc, const long &grad_E_Nmc);

    void stochasticReconfigurationOptimization(const long &Nmc, const double stepSize = 1., const bool flag_dgrad = false); // calc&use gradient error?

    void adamOptimization(const long &Nmc, const bool useSR = false, const bool useGradientError = false, const size_t &max_n_const_values = 20, const bool useAveraging = false,
        const double &lambda = 0, const double &alpha = 0.001, const double &beta1 = 0.9, const double &beta2 = 0.999, const double &epsilon = 10e-8);

    void simulatedAnnealingOptimization(const long &Nmc, const double &iota, const double &kappa, const double &lambda, gsl_siman_params_t &params);

    void nmsimplexOptimization(const long &Nmc, const double &iota, const double &kappa, const double &lambda, const double &rstart = 1.0, const double &rend = 0.01, const size_t &max_n_iter = 0);
};


#endif