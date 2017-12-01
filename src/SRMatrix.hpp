#ifndef SR_MATRIX
#define SR_MATRIX

#include "MCIObservableFunctionInterface.hpp"
#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"

// SR matrix c_ij:
//     c_ij = ⟨Oi Oj⟩
// c_ij is a nxn matrix. We map it to an array according to this rule:
// obs[i*n + j] = c_ij
class SRMatrix: public MCIObservableFunctionInterface
{
protected:
   WaveFunction * _wf;
      
public:
   SRMatrix(WaveFunction * wf):
   MCIObservableFunctionInterface(wf->getNDim(), wf->getNVP() * wf->getNVP()){
      _wf = wf;
   }
   
   virtual ~SRMatrix(){ }


   void observableFunction(const double * in, double *out)
   {
      for (int i=0; i<_wf->getNVP(); ++i){
         for (int j=0; j<_wf->getNVP(); ++j){
            out[i*_wf->getNVP() + j] = _wf->vd1(i,in) * _wf->vd1(j,in);
         }
      }
   }

};


#endif
