#ifndef VMC_STOCHASTICRECONFIGURATIONMCOBSERVABLE_HPP
#define VMC_STOCHASTICRECONFIGURATIONMCOBSERVABLE_HPP

#include "mci/ObservableFunctionInterface.hpp"
#include "vmc/Hamiltonian.hpp"
#include "vmc/WaveFunction.hpp"

namespace vmc
{

class StochasticReconfigurationMCObservable: public mci::ObservableFunctionInterface
{
protected:
    WaveFunction * const _wf;
    Hamiltonian * const _H;

    mci::ObservableFunctionInterface * _clone() const final
    {
        return new StochasticReconfigurationMCObservable(_wf, _H); // ownership needs fix
    }
public:
    StochasticReconfigurationMCObservable(WaveFunction * wf, Hamiltonian * H):
            mci::ObservableFunctionInterface(H->getNDim(), 2*wf->getNVP() + wf->getNVP()*wf->getNVP()),
            _wf(wf), _H(H) {}

    ~StochasticReconfigurationMCObservable() override = default;


    // mci::ObservableFunctionInterface implementation
    void observableFunction(const double * in, const mci::SamplingFunctionContainer &pdfs, double * out) override
    {
        // out is made in this way (nvp is the number of variational parameters):
        // out[0:nvp-1] = Oi
        // out[nvp:2*nvp-1] = HOi
        // out[2*nvp:2*nvp+nvp*nvp-1] = OiOj    according to the rule OiOj[i, j] = OiOj[j + i*nvp]

        // number of variational parameters
        const int nvp = _wf->getNVP();
        // local energy
        const double Hloc = _H->localPBKineticEnergy(in) + _H->localPotentialEnergy(in);

        // store the elements Oi and HOi
        for (int i = 0; i < nvp; ++i) {
            out[i] = _wf->getVD1DivByWF(i);    // Oi
            out[i + nvp] = Hloc*out[i];    //HOi
        }
        // store the elements OiOj
        for (int i = 0; i < nvp; ++i) {
            for (int j = 0; j < nvp; ++j) {
                out[2*nvp + i*nvp + j] = out[i]*out[j];   // OiOj
            }
        }
    }
};
} // namespace vmc

#endif
