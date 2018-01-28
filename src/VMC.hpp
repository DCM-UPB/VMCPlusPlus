#ifndef VMC_CLASS
#define VMC_CLASS

#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "ConjugateGradientOptimization.hpp"
#include "StochasticReconfigurationOptimization.hpp"
#include "SimulatedAnnealingOptimization.hpp"

#include "MCIntegrator.hpp"

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
      if (wf->getNDim() != H->getNDim())
         throw std::invalid_argument( "Error VMC: ndim different between wf and H" );
      _wf=wf; 
      _H=H;
      _mci = new MCI(_H->getNDim());
   }

   ~VMC(){
      delete _mci;
   }
   
   
   // Monte Carlo Integral within VMC should be performed using the MCI object provided by VMC
   MCI * getMCI(){return _mci;}
   
   
   // Computation of the variational energy
   void computeVariationalEnergy(const long & Nmc, double * E, double * dE);
   

   // Wave Function Optimization Methods
   void conjugateGradientOptimization(const long &E_Nmc, const long &grad_E_Nmc);
   
   void stochasticReconfigurationOptimization(const long &Nmc);
   
   void simulatedAnnealingOptimization(const long &Nmc, const double &iota, const double &kappa, const double &lambda, gsl_siman_params_t &params);

};


#endif
