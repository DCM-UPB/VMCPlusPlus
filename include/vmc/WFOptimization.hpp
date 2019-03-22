#ifndef VMC_WFOPTIMIZATION_HPP
#define VMC_WFOPTIMIZATION_HPP

#include "mci/MCIntegrator.hpp"
#include "vmc/Hamiltonian.hpp"
#include "vmc/WaveFunction.hpp"


class WFOptimization{
protected:
    WaveFunction * const _wf;
    Hamiltonian * const _H;
    mci::MCI * const _mci;

public:
    WFOptimization(WaveFunction * wf, Hamiltonian * H, mci::MCI * mci):
        _wf(wf), _H(H), _mci(mci) {}

    virtual ~WFOptimization()= default;

    // getters
    WaveFunction * getWF(){return _wf;}
    Hamiltonian * getH(){return _H;}
    mci::MCI * getMCI(){return _mci;};

    // optimize the wf
    virtual void optimizeWF() = 0;
};



#endif
