#ifndef WF_OPTIMIZATION
#define WF_OPTIMIZATION

#include "vmc/WaveFunction.hpp"
#include "vmc/Hamiltonian.hpp"
#include "mci/MCIntegrator.hpp"


class WFOptimization{
protected:
    WaveFunction * _wf;
    Hamiltonian * _H;
    MCI * _mci;

public:
    WFOptimization(WaveFunction * wf, Hamiltonian * H, MCI * mci){
        _wf=wf;
        _H=H;
        _mci = mci;
    }

    virtual ~WFOptimization(){}

    // getters
    WaveFunction * getWF(){return _wf;}
    Hamiltonian * getH(){return _H;}
    MCI * getMCI(){return _mci;};

    // optimize the wf
    virtual void optimizeWF() = 0;
};



#endif
