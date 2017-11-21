#include "VMC.hpp"

#include "ConjGrad.hpp"


// --- Optimization

void VMC::conjugateGradientOptimization()
{
   // declare the Conjugate Gradient object
   ConjGrad cjgrad(this);
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

void VMC::f(const double * in, double &f, double &df)
{
   //set variational parameters in the wave function
   _wf->setVP(in);

   double * E = new double[4];
   double * dE = new double[4];

   this->computeEnergy(_ENmc, E, dE);

   f=E[0]; df=dE[0];

   //using namespace std;
   //cout << "E = " << E[0] << " +- " << dE[0] << endl;

   delete[] E;
   delete[] dE;
}


void VMC::grad(const double *in, double *g, double *dg)
{
   //set variational parameters in the wave function
   _wf->setVP(in);

   this->computeEnergyGradient(_GNmc, g, dg);
}


// --- Public methods

void VMC::computeEnergy(const long & Nmc, double * E, double * dE)
{
   _mcenergy->clearSamplingFunctions();
   _mcenergy->addSamplingFunction(_wf);
   _mcenergy->clearObservables();
   _mcenergy->addObservable(_H);
   _mcenergy->integrate(Nmc,E,dE);
}


void VMC::computeEnergyGradient(const long & Nmc, double *gradE, double * dgradE)
{
   _mcenergy->clearSamplingFunctions();
   _mcenergy->addSamplingFunction(_wf);
   _mcenergy->clearObservables();
   _mcenergy->addObservable(_H);
   _mcenergy->addObservable(_VG);
   double * obs = new double[4+2*_wf->getNVP()];
   double * dobs = new double[4+2*_wf->getNVP()];
   _mcenergy->integrate(Nmc,obs,dobs);
   for (int i=0; i<_wf->getNVP(); ++i)
   {
      gradE[i] = 2.*( obs[4+_wf->getNVP()+i] - obs[0]*obs[4+i] );
      dgradE[i] = 2.*( dobs[4+_wf->getNVP()+i] + (obs[0]*obs[4+i])*(dobs[0]/obs[0]+dobs[4+i]/obs[4+i]) ) ;
   }
   delete[] obs;
   delete[] dobs;
}

