#ifndef STOCHASTIC_RECONFIGURATION_MC_OBSERVABLE
#define STOCHASTIC_RECONFIGURATION_MC_OBSERVABLE

#include "MCIObservableFunctionInterface.hpp"
#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"



class StochasticReconfigurationMCObservable: public MCIObservableFunctionInterface
{
protected:
    WaveFunction * _wf;
    Hamiltonian * _H;

public:
    StochasticReconfigurationMCObservable(WaveFunction * wf, Hamiltonian * H):
        MCIObservableFunctionInterface(H->getNDim(), 2*wf->getNVP() + wf->getNVP()*wf->getNVP()){
        _wf = wf;
        _H = H;
    }

    virtual ~StochasticReconfigurationMCObservable(){}


    // MCIObservableFunctionInterface implementation
    void observableFunction(const double * in, double * out){
        // out is made in this way (nvp is the number of variational parameters):
        // out[0:nvp-1] = Oi
        // out[nvp:2*nvp-1] = HOi
        // out[2*nvp:2*nvp+nvp*nvp-1] = OiOj    according to the rule OiOj[i, j] = OiOj[j + i*nvp]

        // number of variational parameters
        const int nvp = _wf->getNVP();
        // local energy
        const double Hloc = _H->localPBKineticEnergy(in) + _H->localPotentialEnergy(in);

        // variational derivatives
        double * vd1 = new double[nvp];
        for (int i=0; i<nvp; ++i){
            vd1[i] = _wf->getVD1LogWF(i);
        }
        // store the elements Oi and HOi
        for (int i=0; i<nvp; ++i){
            out[i] = vd1[i];    // Oi
            out[i+nvp] = Hloc * out[i];    //HOi
        }
        // store the elements OiOj
        for (int i=0; i<nvp; ++i){
            for (int j=0; j<nvp; ++j){
                out[2*nvp + j + i*nvp] = vd1[i]*vd1[j];   // OiOj
            }
        }
        // free resources
        delete[] vd1;
    }

};


#endif
