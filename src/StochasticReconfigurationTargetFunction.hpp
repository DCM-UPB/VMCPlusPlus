#ifndef STOCHASTIC_RECONFIGURATION_TARGET_FUNCTION
#define STOCHASTIC_RECONFIGURATION_TARGET_FUNCTION

#include "StochasticReconfigurationMCObservable.hpp"
#include "MCIntegrator.hpp"
#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "NoisyFunction.hpp"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <iostream>



class StochasticReconfigurationTargetFunction: public NoisyFunctionWithGradient
{
   protected:
      WaveFunction * _wf;
      Hamiltonian * _H;
      long _Nmc;
      MCI * _mci;
      
   public:
      StochasticReconfigurationTargetFunction(WaveFunction * wf, Hamiltonian * H, const long & Nmc, MCI * mci): 
      NoisyFunctionWithGradient(wf->getNVP()){
         _wf = wf;
         _H = H;
         _Nmc = Nmc;
         _mci = mci;
      }
      
      virtual ~StochasticReconfigurationTargetFunction(){}
      
      
      // NoisyFunctionWithGradient implementation
      void f(const double *vp, double &f, double &df){
         std::cout << "   f    " << std::endl;
         
         // set the variational parameters given as input
         _wf->setVP(vp);
         // set up the MC integrator
         _mci->clearSamplingFunctions(); _mci->addSamplingFunction(_wf);
         _mci->clearObservables(); _mci->addObservable(_H);
         // perform the integral and store the values
         double * obs = new double[4];
         double * dobs = new double[4];
         std::cout << "vp = " << vp[0] << "   " << vp[1] << std::endl;
         std::cout << "integrate " << _Nmc << "   " << _mci->getNSampF() << "   " << _mci->getNObs() << "   X = " << _mci->getX(0) <<  std::endl;
         _mci->integrate(_Nmc, obs, dobs);
         std::cout << "end integrate" << std::endl;
         f = obs[0];
         df = dobs[0];
         // free resources
         delete dobs;
         delete obs;
      }
      
      void grad(const double *vp, double *grad_E, double *dgrad_E){
         std::cout << "   grad    " << std::endl;
         
         // set the variational parameters given as input
         _wf->setVP(vp);
         // set up the MC integrator
         _mci->clearSamplingFunctions(); _mci->addSamplingFunction(_wf);
         _mci->clearObservables();
         _mci->addObservable(_H);
         StochasticReconfigurationMCObservable * mc_obs = new StochasticReconfigurationMCObservable(_wf, _H);
         _mci->addObservable(mc_obs);
         // perform the integral and store the values
         double * obs = new double[4 + 2*_wf->getNVP() + _wf->getNVP()*_wf->getNVP()];
         double * dobs = new double[4 + 2*_wf->getNVP() + _wf->getNVP()*_wf->getNVP()];
         _mci->integrate(_Nmc, obs, dobs);
         // create pointers for ease of use and readability
         double * H = obs;
         double * dH = dobs;
         double * Oi = obs+4;
         double * dOi = dobs+4;
         double * HOi = obs+4+_wf->getNVP();
         double * dHOi = dobs+4+_wf->getNVP();
         double * OiOj = obs+4+2*_wf->getNVP();
         double * dOiOj = dobs+4+2*_wf->getNVP();
         
         // --- compute direction (or gradient) to follow
         gsl_matrix * sij = gsl_matrix_alloc(_wf->getNVP(), _wf->getNVP());
         gsl_matrix * rdsij = gsl_matrix_alloc(_wf->getNVP(), _wf->getNVP());   // relative error, i.e. error/value
         for (int i=0; i<_wf->getNVP(); ++i){
            for (int j=0; j<_wf->getNVP(); ++j){
               gsl_matrix_set(sij, i, j, OiOj[i*_wf->getNVP() + j] - Oi[i] * Oi[j]);
               gsl_matrix_set(rdsij, i, j, 
                  (dOiOj[i*_wf->getNVP() + j] + abs(Oi[i]*Oi[j])*( (dOi[i]/Oi[i]) + (dOi[j]/Oi[j]) ))  
                     /  gsl_matrix_get(sij, i, j) );
            }
         }
         gsl_vector * fi = gsl_vector_alloc(_wf->getNVP());
         gsl_vector * rdfi = gsl_vector_alloc(_wf->getNVP());   // relative error, i.e. error/value
         for (int i=0; i<_wf->getNVP(); ++i){
            gsl_vector_set(fi, i, H[0]*Oi[i] - HOi[i]);
            gsl_vector_set(rdfi, i, 
               (abs(H[0]*Oi[i])*( (dH[0]/H[0]) + (dOi[i]/Oi[i]) ) + dHOi[i])
                  / gsl_vector_get(fi, i)  );
         }
         // invert matrix using SVD
         const double SVD_MIN = 1.0e-9;
         // matrix and vectors needed for the SVD
         gsl_matrix * V = gsl_matrix_alloc(_wf->getNVP(), _wf->getNVP());
         gsl_vector * S = gsl_vector_alloc(_wf->getNVP());
         gsl_vector * work = gsl_vector_alloc(_wf->getNVP());
         // run the Single Value Decomposition
         gsl_linalg_SV_decomp(sij, V, S, work);
         // assemble the inverse matrix
         gsl_matrix * Isij = gsl_matrix_alloc(_wf->getNVP(), _wf->getNVP());
         for (int i=0; i<_wf->getNVP(); ++i){
            for (int j=0; j<_wf->getNVP(); ++j){
               gsl_matrix_set(Isij, i, j, 0.);
            }
         }
         for (int i=0; i<_wf->getNVP(); ++i){
            if (gsl_vector_get(S, i) > SVD_MIN * gsl_vector_get(S, 0)){
               gsl_matrix_set(Isij, i, i,   1./gsl_vector_get(S, i)  );
            } else {
               gsl_matrix_set(Isij, i, i,   0.  );
            }
         }
         gsl_matrix * mm = gsl_matrix_alloc(_wf->getNVP(), _wf->getNVP());
         gsl_matrix_transpose(V);
         gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Isij, V, 0.0, mm);
         gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, sij, mm, 0.0, Isij);         
         // --- finally find the direction to follow
         double foo;
         for (int i=0; i<_wf->getNVP(); ++i){
            grad_E[i] = 0.;
            dgrad_E[i] = 0.;
            for (int k=0; k<_wf->getNVP(); ++k){
               foo = gsl_vector_get(fi, k)*gsl_matrix_get(Isij, k, i);
               grad_E[i] -= foo;  // the minus sign is there beacuse the NoisyFunMin library will follow 'minus the gradient'
               dgrad_E[i] += abs(foo) * ( gsl_vector_get(rdfi, k) + gsl_matrix_get(rdsij, k, i) );  // not correct, just a rough estimation
            }
         }
      
         // free resources
         gsl_matrix_free(mm);
         gsl_vector_free(work);
         gsl_vector_free(S);
         gsl_matrix_free(V);
         gsl_matrix_free(Isij);
         gsl_vector_free(rdfi);
         gsl_vector_free(fi);
         gsl_matrix_free(rdsij);
         gsl_matrix_free(sij); 
         
         delete[] dobs;
         delete[] obs;
      }

};


#endif
