#include "VMC.hpp"

#include "ConjGrad.hpp"

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>


// --- Optimization

void VMC::conjugateGradientOptimization(){
   // declare the Conjugate Gradient object
   ConjGrad cjgrad(this);
   if (_grad_type == "SR")
      cjgrad.configureToFollowSimpleGradient();
   // allocate an array that will contain the wave function variational parameters
   double * wfpar = new double[_wf->getNVP()];
   // get the variational parameters
   _wf->getVP(wfpar);
   // set the actual variational parameters as starting point for the Conjugate Gradient algorithm
   cjgrad.setX(wfpar);
   // find the optimal parameters by minimizing the energy with the Conjugate Gradient algorithm
   cjgrad.findMin();
   // set the found parameters in the wave function
   cjgrad.getX(wfpar);
   _wf->setVP(wfpar);
   // free memory
   delete[] wfpar;
}



// --- Implementation of the NoisyFunctionWithGradient interface

void VMC::f(const double * in, double &f, double &df){
   //set variational parameters in the wave function
   _wf->setVP(in);

   double * E = new double[4];
   double * dE = new double[4];

   this->computeEnergy(_ENmc, E, dE);

   f=E[0]; df=dE[0];

   delete[] E;
   delete[] dE;
}


void VMC::grad(const double *in, double *g, double *dg){
   //set variational parameters in the wave function
   _wf->setVP(in);

   this->computeEnergyGradient(_GNmc, g, dg);
}


// --- Public methods

void VMC::computeEnergy(const long & Nmc, double * E, double * dE){
   _mcenergy->clearSamplingFunctions();
   _mcenergy->addSamplingFunction(_wf);
   _mcenergy->clearObservables();
   _mcenergy->addObservable(_H);
   _mcenergy->integrate(Nmc,E,dE);
}


void VMC::computeEnergyGradient(const long & Nmc, double *gradE, double * dgradE){
   _mcenergy->clearSamplingFunctions();
   _mcenergy->addSamplingFunction(_wf);
   _mcenergy->clearObservables();
   
   _mcenergy->addObservable(_H);
   int nobs = 4;
   if (_grad_type == "gradE"){
      _mcenergy->addObservable(_VG);
      nobs += 2*_wf->getNVP();
   } else if (_grad_type == "SR"){
      _mcenergy->addObservable(_VG);
      nobs += 2*_wf->getNVP();
      _mcenergy->addObservable(_SRM);
      nobs += _wf->getNVP() * _wf->getNVP();
   }
   
   double * obs = new double[nobs];
   double * dobs = new double[nobs];
   
   _mcenergy->integrate(Nmc,obs,dobs);
   
   if (_grad_type == "gradE"){
      // create pointers for ease of use and readability
      double * H = obs;
      double * dH = dobs;
      double * Oi = obs+4;
      double * dOi = dobs+4;
      double * HOi = obs+4+_wf->getNVP();
      double * dHOi = dobs+4+_wf->getNVP();
      // compute direction (or gradient) to follow
      for (int i=0; i<_wf->getNVP(); ++i){
         gradE[i] = 2.*( HOi[i] - H[0]*Oi[i] );
         dgradE[i] = 2.*( dHOi[i] + abs(H[0]*Oi[i])*(dH[0]/H[0]+dOi[i]/Oi[i]) ) ;
      }
   } else if (_grad_type == "SR"){
      // create pointers for ease of use and readability
      double * H = obs;
      double * dH = dobs;
      double * Oi = obs+4;
      double * dOi = dobs+4;
      double * HOi = obs+4+_wf->getNVP();
      double * dHOi = dobs+4+_wf->getNVP();
      double * OiOj = obs+4+_wf->getNVP()+_wf->getNVP();
      double * dOiOj = dobs+4+_wf->getNVP()+_wf->getNVP();
      // construct sij and fi
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
      // invert matrix sij -> Isij
      gsl_permutation * p = gsl_permutation_alloc(_wf->getNVP());;
      int * signum = new int;
      gsl_linalg_LU_decomp(sij, p, signum);
      gsl_matrix * Isij = gsl_matrix_alloc(_wf->getNVP(), _wf->getNVP());
      gsl_linalg_LU_invert(sij, p, Isij);
      // compute direction (or gradient) to follow
      double foo;
      for (int i=0; i<_wf->getNVP(); ++i){
         gradE[i] = 0.;
         dgradE[i] = 0.;
         for (int k=0; k<_wf->getNVP(); ++k){
            foo = gsl_vector_get(fi, k)*gsl_matrix_get(Isij, k, i);
            gradE[i] -= foo;  // the minus sign is beacuse the NoisyFunMin library will follow 'minus the gradient'
            dgradE[i] += abs(foo) * ( gsl_vector_get(rdfi, k) + gsl_matrix_get(rdsij, k, i) );
         }
      }
      // free memory
      gsl_matrix_free(Isij);
      delete signum;
      gsl_permutation_free(p);
      gsl_vector_free(rdfi);
      gsl_vector_free(fi);
      gsl_matrix_free(rdsij);
      gsl_matrix_free(sij);
   }
   
   
   // debugging code
   //using namespace std;
   //cout << "gradE = " << gradE[0] << " +- " << dgradE[0] << "    " << gradE[1] << " +- " << dgradE[1] << endl;
   //cout << "    ||gradE|| = " << sqrt(gradE[0]*gradE[0]+gradE[1]*gradE[1]) << " +- " << sqrt(gradE[0]*gradE[0]+gradE[1]*gradE[1]) << endl << endl;
   
   
   delete[] obs;
   delete[] dobs;
}

