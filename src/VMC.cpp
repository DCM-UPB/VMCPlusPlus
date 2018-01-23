#include "VMC.hpp"

#include "MCIntegrator.hpp"



void VMC::computeVariationalEnergy(const long & Nmc, double * E, double * dE){
   getMCI()->clearSamplingFunctions(); getMCI()->addSamplingFunction(_wf);
   getMCI()->clearObservables(); getMCI()->addObservable(_H);
   getMCI()->integrate(Nmc, E, dE);
}



//#include "ConjGrad.hpp"
//#include "DynamicDescent.hpp"
//
//#include <sstream>
//#include <math.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_permutation.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_blas.h>
//
//
//
//// --- WF Optimization
//
//
//
//
//
//
//
//// GSL Simulated Annelaing
//
//int _NVP;   // number of dimensions
//double _minE;
//double * _optP;
//VMC * _simanVMC;   // VMC static object, needed for computing the enrgy
//
//double _simanEnergy(void *xp){
//   using namespace std;
//   // extract variational parameters
//   double * p = ((double *) xp);
//   
//   // compute energy
//   double E, dE;
//   _simanVMC->f(p, E, dE);
//   //cout << "   ---> E = " << E << "(" << dE << ")";
//   
//   // compute normalization
//   double norm = 0.;
//   for (int i=0; i<_NVP; ++i){
//      norm += p[i]*p[i];
//   }
//   //cout << "    norm = " << norm << endl;
//   
//   // return
//   return dE; // + 1.0 * norm;
//}
//
//
//void _simanStep(const gsl_rng * r, void * xp, double step_size){
//   // extract old variational parameters
//   double * p = ((double *) xp);
//   // compute new variational parameters
//   double eta;
//   const double double_step_size = 2. * step_size;
//   for (int i=0; i<_NVP; ++i){
//      eta = gsl_rng_uniform(r);
//      p[i] = p[i] + (eta-0.5) * double_step_size;
//   }
//}
//
//
//double _simanDistance(void *xp, void *yp){
//   // extract the two variatonal parameters of which I must compute the distance
//   double * x = ((double *) xp);
//   double * y = ((double *) yp);
//   // compute the distance
//   double dist = 0.;
//   for (int i=0; i<_NVP; ++i){
//      dist += (y[i] - x[i]) * (y[i] - x[i]);
//   }
//   dist = sqrt(dist);
//   // return
//   return dist;
//}
//
//
//void _simanPrint(void *xp){
//   using namespace std;
//   // extract variational parameters
//   double * p = ((double *) xp);
//   // print current parameters
//   cout << endl << "E = " << _minE << "     p ---> ";
//   for (int i=0; i<_NVP; ++i){
//      cout << p[i] << "    ";
//   }
//   cout << endl;
//}
//
//
//
//// --- Optimization
//
//void VMC::conjugateGradientOptimization(){
//   // declare the Conjugate Gradient object
//   ConjGrad cjgrad(this);
//   if (_grad_type == "SR")
//      cjgrad.configureToFollowSimpleGradient();
//   // allocate an array that will contain the wave function variational parameters
//   double * wfpar = new double[_wf->getNVP()];
//   // get the variational parameters
//   _wf->getVP(wfpar);
//   // set the actual variational parameters as starting point for the Conjugate Gradient algorithm
//   cjgrad.setX(wfpar);
//   // find the optimal parameters by minimizing the energy with the Conjugate Gradient algorithm
//   cjgrad.findMin();
//   // set the found parameters in the wave function
//   cjgrad.getX(wfpar);
//   _wf->setVP(wfpar);
//   // free memory
//   delete[] wfpar;
//}
//
//
//void VMC::dynamicDescentOptimization(){
//   // declare the Dynamic Descent object
//   DynamicDescent dynamdesc(this);
//   // allocate an array that will contain the wave function variational parameters
//   double * wfpar = new double[_wf->getNVP()];
//   // get the variational parameters
//   _wf->getVP(wfpar);
//   // set the actual variational parameters as starting point for the Conjugate Gradient algorithm
//   dynamdesc.setX(wfpar);
//   // find the optimal parameters by minimizing the energy with the Conjugate Gradient algorithm
//   dynamdesc.findMin();
//   // set the found parameters in the wave function
//   dynamdesc.getX(wfpar);
//   _wf->setVP(wfpar);
//   // free memory
//   delete[] wfpar;
//}
//
//
//
//void VMC::simulatedAnnealingOptimization(){
//   using namespace std;
//   
//   // simulated annealing constants
//   const int N_TRIES = 7;
//   const int ITERS_FIXED_T = 7;
//   const double STEP_SIZE = 1.;
//   const double K = 1.;
//   const double T_INITIAL = 100.;
//   const double MU_T = 1.01;
//   const double T_MIN = 0.01;
//   gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};
//   
//   // random generator used by the simulated annealing 
//   const gsl_rng_type * T;
//   T = gsl_rng_default;
//   gsl_rng * r;
//   r = gsl_rng_alloc(T);
//   
//   // get the initial parameters
//   double * wfpar = new double[_wf->getNVP()];
//   _wf->getVP(wfpar);
//   
//   _NVP = _wf->getNVP();
//   _optP = new double[_NVP];
//   _minE = 1000.;
//   _simanVMC = this;
//   gsl_siman_solve(r, wfpar, _simanEnergy, _simanStep, _simanDistance, _simanPrint, NULL, NULL, NULL, _wf->getNVP()*sizeof(double), params);
//   
//   // // store the new parameters 
//   // _wf->setVP(_optP);
//   // double E, dE;
//   // this->computeEnergy(_ENmc, &E, &dE);
//   // cout << "newly estimated energy = " << E << " (" << dE << ")" << endl;
//   
//   // free resources
//   delete[] _optP;
//   delete[] wfpar;
//   gsl_rng_free(r);
//}
//
//
//
//
//// --- Implementation of the NoisyFunctionWithGradient interface
//
//void VMC::f(const double * in, double &f, double &df){
//   //set variational parameters in the wave function
//   _wf->setVP(in);
//
//   double * E = new double[4];
//   double * dE = new double[4];
//
//   this->computeEnergy(_ENmc, E, dE);
//
//   f=E[0]; df=dE[0];
//
//   delete[] E;
//   delete[] dE;
//}
//
//
//void VMC::grad(const double *in, double *g, double *dg){
//   //set variational parameters in the wave function
//   _wf->setVP(in);
//
//   this->computeEnergyGradient(_GNmc, g, dg);
//}
//
//
//// --- Public methods
//
//void VMC::computeEnergy(const long & Nmc, double * E, double * dE){
//   _mcenergy->clearSamplingFunctions();
//   _mcenergy->addSamplingFunction(_wf);
//   _mcenergy->clearObservables();
//   _mcenergy->addObservable(_H);
//   _mcenergy->integrate(Nmc,E,dE);
//}
//
//
//void VMC::computeEnergyGradient(const long & Nmc, double *gradE, double * dgradE){
//   using namespace std;
//   
//   _mcenergy->clearSamplingFunctions();
//   _mcenergy->addSamplingFunction(_wf);
//   _mcenergy->clearObservables();
//   
//   _mcenergy->addObservable(_H);
//   int nobs = 4;
//   if (_grad_type == "gradE"){
//      _mcenergy->addObservable(_VG);
//      nobs += 2*_wf->getNVP();
//   } else if (_grad_type == "SR"){
//      _mcenergy->addObservable(_VG);
//      nobs += 2*_wf->getNVP();
//      _mcenergy->addObservable(_SRM);
//      nobs += _wf->getNVP() * _wf->getNVP();
//   }
//   
//   double * obs = new double[nobs];
//   double * dobs = new double[nobs];
//   
//   cout << "status 1" << endl;
//   _mcenergy->integrate(Nmc,obs,dobs);
//   cout << "status 2" << endl;
//   
//   if (_grad_type == "gradE"){
//      // create pointers for ease of use and readability
//      double * H = obs;
//      double * dH = dobs;
//      double * Oi = obs+4;
//      double * dOi = dobs+4;
//      double * HOi = obs+4+_wf->getNVP();
//      double * dHOi = dobs+4+_wf->getNVP();
//      // compute direction (or gradient) to follow
//      for (int i=0; i<_wf->getNVP(); ++i){
//         gradE[i] = 2.*( HOi[i] - H[0]*Oi[i] );
//         dgradE[i] = 2.*( dHOi[i] + abs(H[0]*Oi[i])*(dH[0]/H[0]+dOi[i]/Oi[i]) ) ;
//      }
//   } else if (_grad_type == "SR"){
//      // create pointers for ease of use and readability
//      double * H = obs;
//      double * dH = dobs;
//      double * Oi = obs+4;
//      double * dOi = dobs+4;
//      double * HOi = obs+4+_wf->getNVP();
//      double * dHOi = dobs+4+_wf->getNVP();
//      double * OiOj = obs+4+_wf->getNVP()+_wf->getNVP();
//      double * dOiOj = dobs+4+_wf->getNVP()+_wf->getNVP();
//      // construct sij and fi
//      gsl_matrix * sij = gsl_matrix_alloc(_wf->getNVP(), _wf->getNVP());
//      gsl_matrix * rdsij = gsl_matrix_alloc(_wf->getNVP(), _wf->getNVP());   // relative error, i.e. error/value
//      for (int i=0; i<_wf->getNVP(); ++i){
//         for (int j=0; j<_wf->getNVP(); ++j){
//            gsl_matrix_set(sij, i, j, OiOj[i*_wf->getNVP() + j] - Oi[i] * Oi[j]);
//            gsl_matrix_set(rdsij, i, j, 
//               (dOiOj[i*_wf->getNVP() + j] + abs(Oi[i]*Oi[j])*( (dOi[i]/Oi[i]) + (dOi[j]/Oi[j]) ))  
//                  /  gsl_matrix_get(sij, i, j) );
//         }
//      }
//      gsl_vector * fi = gsl_vector_alloc(_wf->getNVP());
//      gsl_vector * rdfi = gsl_vector_alloc(_wf->getNVP());   // relative error, i.e. error/value
//      for (int i=0; i<_wf->getNVP(); ++i){
//         gsl_vector_set(fi, i, H[0]*Oi[i] - HOi[i]);
//         gsl_vector_set(rdfi, i, 
//            (abs(H[0]*Oi[i])*( (dH[0]/H[0]) + (dOi[i]/Oi[i]) ) + dHOi[i])
//               / gsl_vector_get(fi, i)  );
//      }
//      
//      
//      
//      // --- INVERT MATRIX - OLD
//      // // invert matrix sij -> Isij
//      // gsl_permutation * p = gsl_permutation_alloc(_wf->getNVP());;
//      // int * signum = new int;
//      // gsl_linalg_LU_decomp(sij, p, signum);
//      // gsl_matrix * Isij = gsl_matrix_alloc(_wf->getNVP(), _wf->getNVP());
//      // gsl_linalg_LU_invert(sij, p, Isij);
//      //
//      //delete signum;
//      //gsl_permutation_free(p);
//      
//      
//      
//      // --- INVERT MATRIX USING SVD
//      const double SVD_MIN = 1.0e-9;
//      // matrix and vectors needed for the SVD
//      gsl_matrix * V = gsl_matrix_alloc(_wf->getNVP(), _wf->getNVP());
//      gsl_vector * S = gsl_vector_alloc(_wf->getNVP());
//      gsl_vector * work = gsl_vector_alloc(_wf->getNVP());
//      // run the Single Value Decomposition
//      gsl_linalg_SV_decomp(sij, V, S, work);
//      // assemble the inverse matrix
//      gsl_matrix * Isij = gsl_matrix_alloc(_wf->getNVP(), _wf->getNVP());
//      for (int i=0; i<_wf->getNVP(); ++i){
//         for (int j=0; j<_wf->getNVP(); ++j){
//            gsl_matrix_set(Isij, i, j, 0.);
//         }
//      }
//      for (int i=0; i<_wf->getNVP(); ++i){
//         if (gsl_vector_get(S, i) > SVD_MIN * gsl_vector_get(S, 0)){
//            gsl_matrix_set(Isij, i, i,   1./gsl_vector_get(S, i)  );
//         } else {
//            gsl_matrix_set(Isij, i, i,   0.  );
//         }
//      }
//      gsl_matrix * mm = gsl_matrix_alloc(_wf->getNVP(), _wf->getNVP());
//      gsl_matrix_transpose(V);
//      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Isij, V, 0.0, mm);
//      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, sij, mm, 0.0, Isij);
//      // free temporary memory
//      gsl_matrix_free(mm);
//      gsl_vector_free(work);
//      gsl_vector_free(S);
//      gsl_matrix_free(V);
//      
//      
//      // compute direction (or gradient) to follow
//      double foo;
//      for (int i=0; i<_wf->getNVP(); ++i){
//         gradE[i] = 0.;
//         dgradE[i] = 0.;
//         for (int k=0; k<_wf->getNVP(); ++k){
//            foo = gsl_vector_get(fi, k)*gsl_matrix_get(Isij, k, i);
//            gradE[i] -= foo;  // the minus sign is beacuse the NoisyFunMin library will follow 'minus the gradient'
//            dgradE[i] += abs(foo) * ( gsl_vector_get(rdfi, k) + gsl_matrix_get(rdsij, k, i) );
//         }
//      }
//      
//      // free memory
//      gsl_matrix_free(Isij);
//      gsl_vector_free(rdfi);
//      gsl_vector_free(fi);
//      gsl_matrix_free(rdsij);
//      gsl_matrix_free(sij);      
//   }
//   
//   
//   // debugging code
//   //using namespace std;
//   //cout << "gradE = " << gradE[0] << " +- " << dgradE[0] << "    " << gradE[1] << " +- " << dgradE[1] << endl;
//   //cout << "    ||gradE|| = " << sqrt(gradE[0]*gradE[0]+gradE[1]*gradE[1]) << " +- " << sqrt(gradE[0]*gradE[0]+gradE[1]*gradE[1]) << endl << endl;
//   
//   
//   delete[] obs;
//   delete[] dobs;
//   
//}

