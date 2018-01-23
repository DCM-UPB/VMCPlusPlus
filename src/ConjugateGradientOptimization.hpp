#ifndef CONJUGATE_GRADIENT_OPTIMIZATION
#define CONJUGATE_GRADIENT_OPTIMIZATION

#include "WFOptimization.hpp"
#include "ConjugateGradientTargetFunction.hpp"
#include "ConjGrad.hpp"

#include "MCIntegrator.hpp"


class ConjugateGradientOptimization: public WFOptimization{
   
private:
   long _E_Nmc;
   long _grad_E_Nmc;
   MCI * _mci;

public:
   ConjugateGradientOptimization(WaveFunction * wf, Hamiltonian * H, const long &E_Nmc, const long &grad_E_Nmc, MCI * mci): WFOptimization(wf, H, mci){
      _E_Nmc = E_Nmc; _grad_E_Nmc = grad_E_Nmc;
   }
   virtual ~ConjugateGradientOptimization(){}
   
   // optimization
   void optimizeWF(){
      // create targetfunction
      ConjugateGradientTargetFunction * targetf = new ConjugateGradientTargetFunction(_wf, _H, _E_Nmc, _grad_E_Nmc, getMCI());
      // declare the Conjugate Gradient object
      ConjGrad * cjgrad = new ConjGrad(targetf);
      // allocate an array that will contain the wave function variational parameters
      double * wfpar = new double[_wf->getNVP()];
      // get the variational parameters
      _wf->getVP(wfpar);
      // set the actual variational parameters as starting point for the Conjugate Gradient algorithm
      cjgrad->setX(wfpar);
      // find the optimal parameters by minimizing the energy with the Conjugate Gradient algorithm
      cjgrad->findMin();
      // set the found parameters in the wave function
      cjgrad->getX(wfpar);
      _wf->setVP(wfpar);
      // free memory
      delete[] wfpar;
      delete cjgrad;
      delete targetf;
   }
   
};



#endif