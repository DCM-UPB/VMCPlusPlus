#ifndef SIMULATED_ANNEALING_OPTIMIZATION
#define SIMULATED_ANNEALING_OPTIMIZATION

#include "WFOptimization.hpp"

#include <gsl/gsl_siman.h>
#include <cmath>



namespace vmc_siman{
   
   WaveFunction * wf;
   Hamiltonian * H;
   long E_Nmc;
   MCI * mci;
   
   // Simulated Annealing target function parameters
   double iota;   // factor that applies to the energy, for the target function
   double kappa;   // factor that applies to the energy's standard deviation, for the target function
   double lambda;  // normalization factor, for the target function


   double simanTarget(void *xp){
      using namespace std;
   
      // extract the variational parameters and apply them to the wf
      double * x = ((double *) xp);
      // apply the parameters to the wf
      wf->setVP(x);
      
      // compute the energy and its standard deviation
      double * energy = new double[4]; // energy
      double * d_energy = new double[4]; // energy error bar
      mci->clearSamplingFunctions(); mci->addSamplingFunction(wf);
      mci->clearObservables(); mci->addObservable(H);
      mci->integrate(E_Nmc, energy, d_energy);
   
      // compute the normalization
      double norm = 0.;
      for (int i=0; i<wf->getNVP(); ++i)
         norm += x[i]*x[i];
      norm = sqrt(norm)/wf->getNVP();
      
      // assemble the target function
      const double target = iota * energy[0] + kappa * d_energy[0] + lambda * norm;
      
      // free resources
      delete[] d_energy;
      delete[] energy;
   
      return target;
   }
   
   void simanStep(const gsl_rng * r, void * xp, double step_size){   
      // extract old variational parameters
      double * x = ((double *) xp);
      
      // compute new variational parameters
      double eta;
      const double double_step_size = 2. * step_size;
      for (int i=0; i<wf->getNVP(); ++i){
         eta = gsl_rng_uniform(r);
         x[i] = x[i] + (eta-0.5) * double_step_size;
      }
   }
   
   double simanDistance(void *xp, void *yp){
      // extract the two variatonal parameters of which I must compute the distance
      double * x = ((double *) xp);
      double * y = ((double *) yp);
      
      // compute the distance
      double dist = 0.;
      for (int i=0; i<wf->getNVP(); ++i){
         dist += (y[i] - x[i]) * (y[i] - x[i]);
      }
      dist = sqrt(dist);

      return dist;
   }

}



class SimulatedAnnealingOptimization: public WFOptimization{
private:
   // simulated annealing parameters
   //int N_TRIES = 20;
   //int ITERS_FIXED_T = 20;
   //double STEP_SIZE = 1.;
   //double K = 1.;
   //double T_INITIAL = 20.;
   //double MU_T = 1.1;
   //double T_MIN = 0.01;
   gsl_siman_params_t _params; // = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};

public:
   SimulatedAnnealingOptimization(WaveFunction * wf, Hamiltonian * H, const long &E_Nmc, MCI * mci, const double &iota, const double &kappa, const double &lambda, gsl_siman_params_t &params): WFOptimization(wf, H, mci){
      vmc_siman::wf = wf;
      vmc_siman::H = H;
      vmc_siman::E_Nmc = E_Nmc;
      vmc_siman::mci = mci;
      vmc_siman::iota = iota;
      vmc_siman::kappa = kappa;
      vmc_siman::lambda = lambda;
      _params = params;
   }
   virtual ~SimulatedAnnealingOptimization(){}
   
   // optimization
   void optimizeWF(){  
      // random generator used by the simulated annealing 
      gsl_rng_env_setup();
      const gsl_rng_type * T;
      T = gsl_rng_default;
      gsl_rng * r;
      r = gsl_rng_alloc(T);
      
      // set the initial parameters
      double * vp = new double[_wf->getNVP()];
      _wf->getVP(vp);
      
      // run the simulated annealing solver
      gsl_siman_solve(r, vp, vmc_siman::simanTarget, vmc_siman::simanStep, vmc_siman::simanDistance, NULL, NULL, NULL, NULL, _wf->getNVP()*sizeof(double), _params);
      
      // set the optimal variational parameters
      _wf->setVP(vp);
      
      // free resources
      delete[] vp;
      gsl_rng_free(r);
   }
   
};



#endif