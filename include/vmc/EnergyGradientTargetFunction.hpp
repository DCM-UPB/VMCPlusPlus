#ifndef VMC_ENERGYGRADIENTTARGETFUNCTION_HPP
#define VMC_ENERGYGRADIENTTARGETFUNCTION_HPP


#include "vmc/Hamiltonian.hpp"
#include "vmc/WaveFunction.hpp"

#include "mci/MCIntegrator.hpp"
#include "nfm/NoisyFunction.hpp"

class EnergyGradientTargetFunction: public NoisyFunctionWithGradient
{
protected:
    WaveFunction * const _wf;
    Hamiltonian * const _H;
    mci::MCI * const _mci;
    const int _E_Nmc;
    const int _grad_E_Nmc;
    const double _lambda_reg;

public:
    EnergyGradientTargetFunction(WaveFunction * wf, Hamiltonian * H, const int & E_Nmc, const int &grad_E_Nmc, mci::MCI * mci, const double &lambda_reg = 0.):
        NoisyFunctionWithGradient(wf->getNVP()),
        _wf(wf), _H(H), _mci(mci), _E_Nmc(E_Nmc), _grad_E_Nmc(grad_E_Nmc), _lambda_reg(lambda_reg) {}

    ~EnergyGradientTargetFunction() override= default;


    // NoisyFunctionWithGradient implementation
    void f(const double *vp, double &f, double &df) override;
    void grad(const double *vp, double *grad_E, double *dgrad_E) override;
    void fgrad(const double *vp, double &f, double &df, double *grad_E, double *dgrad_E) override;
};


#endif
