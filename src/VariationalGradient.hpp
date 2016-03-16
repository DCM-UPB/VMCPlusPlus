#ifndef VARIATIONAL_GRADIENT
#define VARIATIONAL_GRADIENT

#include "MCIObservableFunctionInterface.hpp"
#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"

// VariationalGradient: obs[0 - wf->getNVP()] = Variational Derivative of the Wave Function   obs[ wf->getNVP() - 2*wf->getNVP()] = Local Energy times the the Variational Derivative of the Wave Function
class VariationalGradient: public MCIObservableFunctionInterface
{
   protected:
      WaveFunction * _wf;
      Hamiltonian * _H;
   public:
      VariationalGradient(WaveFunction * wf, Hamiltonian * H): MCIObservableFunctionInterface(H->getNDim(),2*wf->getNVP())
      {
         _wf = wf;
         _H = H;
      }
      virtual ~VariationalGradient(){ }

      void observableFunction(const double * in, double *out)
      {
         double Hloc;
         Hloc = _H->localPBKineticEnergy(in) + _H->localPotentialEnergy(in);
         for (int i=0; i<_wf->getNVP(); ++i)
         {
            out[i] = _wf->vd1(i,in);
            out[i+_wf->getNVP()] = out[i] * Hloc;
         }
      }

};


#endif
