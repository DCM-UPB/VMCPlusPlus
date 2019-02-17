#ifndef ENERGY_GRADIENT_MC_OBSERVABLE
#define ENERGY_GRADIENT_MC_OBSERVABLE

#include "mci/MCIObservableFunctionInterface.hpp"
#include "vmc/WaveFunction.hpp"
#include "vmc/Hamiltonian.hpp"

#include <iostream>



class EnergyGradientMCObservable: public MCIObservableFunctionInterface
{
protected:
    WaveFunction * _wf;
    Hamiltonian * _H;

public:
    EnergyGradientMCObservable(WaveFunction * wf, Hamiltonian * H):
        MCIObservableFunctionInterface(H->getTotalNDim(), 2*wf->getNVP()){
        _wf = wf;
        _H = H;
    }

    virtual ~EnergyGradientMCObservable(){}



    // MCIObservableFunctionInterface implementation
    void observableFunction(const double * in, double * out){
        // obs[0 : wf->getNVP()-1] = Variational Derivative of the Wave Function
        // obs[ wf->getNVP() : 2*wf->getNVP()-1] = Local Energy times the Variational Derivative of the Wave Function
        double Hloc;

        Hloc = _H->localPBKineticEnergy(in) + _H->localPotentialEnergy(in);
        for (int i=0; i<_wf->getNVP(); ++i){
            out[i] = _wf->getVD1DivByWF(i);
            out[i+_wf->getNVP()] = out[i] * Hloc;
        }
    }

};


#endif
