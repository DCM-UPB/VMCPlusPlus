#ifndef VMC_CLASS
#define VMC_CLASS

#include <iostream>

#include "MCIntegrator.hpp"
#include "NoisyFunction.hpp"

#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "VariationalGradient.hpp"


class VMC: public NoisyFunctionWithGradient
{
   protected:
      long _ENmc;  // number of MC steps used to compute the energy
      long _GNmc;  // number of MC steps used to compute the energy gradient 
      MCI * _mcenergy;
      WaveFunction * _wf;
      Hamiltonian * _H;
      VariationalGradient * _VG;

   public:
      VMC(WaveFunction * wf, Hamiltonian * H, const int &ENmc, const int &GNmc): NoisyFunctionWithGradient(wf->getNVP())
      {
         using namespace std;
      
         if (wf->getNDim() != H->getNDim())
         {
            cout << "ERROR VMC::VMC() : ndim different netween wf and H" << endl;
            exit(1);
         }
         _wf=wf; 
         _H=H;
         _ENmc=ENmc;
         _GNmc=GNmc;
         // MCI to compute the energy
         _mcenergy = new MCI(_H->getNDim());
         // Variational Gradient
         _VG = new VariationalGradient(_wf,_H);
      }

      ~VMC()
      {
         delete _mcenergy;
         delete _VG;
      }

      // --- Getters
      MCI * getEnergyMCI(){return _mcenergy;}

      // Conjugate Gradient minimization
      void conjugateGradientOptimization();

      // Computation of the energy
      void computeEnergy(const long & Nmc, double * E, double * dE);

      // Computation of the energy gradient
      void computeEnergyGradient(const long & Nmc, double * gradE, double * dgradE);

      // Implementation of the NoisyFunctionWithGradient interface
      void f(const double * in, double &f, double &df);
      void grad(const double *in, double *g, double *dg);

};


#endif
