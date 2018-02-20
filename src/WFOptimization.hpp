#ifndef WF_OPTIMIZATION
#define WF_OPTIMIZATION

#include "MCIntegrator.hpp"


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

    ~WFOptimization(){}

    // getters
    MCI * getMCI(){return _mci;};

    // optimize the wf
    virtual void optimizeWF() = 0;

};



#endif
