#ifndef VMC_ENERGYGRADIENTMCOBSERVABLE_HPP
#define VMC_ENERGYGRADIENTMCOBSERVABLE_HPP

#include "mci/ObservableFunctionInterface.hpp"
#include "vmc/Hamiltonian.hpp"
#include "vmc/WaveFunction.hpp"

#include <iostream>



class EnergyGradientMCObservable: public mci::ObservableFunctionInterface
{
protected:
    WaveFunction * const _wf;
    Hamiltonian * const _H;

    mci::ObservableFunctionInterface * _clone() const final {
        return new EnergyGradientMCObservable(_wf, _H);
    }
public:
    EnergyGradientMCObservable(WaveFunction * wf, Hamiltonian * H):
        mci::ObservableFunctionInterface(H->getTotalNDim(), 2*wf->getNVP()), _wf(wf), _H(H) {}

    ~EnergyGradientMCObservable() override= default;


    // mci::ObservableFunctionInterface implementation
    void observableFunction(const double * in, double * out) override{
        // obs[0 : wf->getNVP()-1] = Variational Derivative of the Wave Function
        // obs[ wf->getNVP() : 2*wf->getNVP()-1] = Local Energy times the Variational Derivative of the Wave Function

        const double Hloc = _H->localPBKineticEnergy(in) + _H->localPotentialEnergy(in);
        for (int i=0; i<_wf->getNVP(); ++i){
            out[i] = _wf->getVD1DivByWF(i);
            out[i+_wf->getNVP()] = out[i] * Hloc;
        }
    }
};


#endif
