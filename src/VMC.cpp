#include "vmc/VMC.hpp"
#include "vmc/MPIVMC.hpp"

namespace vmc
{


// --- Constructors

VMC::VMC(std::unique_ptr<WaveFunction> wf, std::unique_ptr<Hamiltonian> H, const int nskip_eg, const int blksize_eg):
        _mci(H->getTotalNDim()), _wf(wf.get()/*remains valid (until destruct)*/), _H(H.get()),
        _nskip_eg(nskip_eg), _blksize_eg(blksize_eg)
{
    if (_wf->getTotalNDim() != _H->getTotalNDim()) {
        throw std::invalid_argument("[VMC] ndim different between wf and H");
    }
    _mci.addSamplingFunction(std::move(wf));
    _mci.addObservable(std::move(H), _blksize_eg, _nskip_eg);
}

VMC::VMC(const WaveFunction &wf, const Hamiltonian &H, const int nskip_eg, const int blksize_eg):
        _mci(H.getTotalNDim()),
        _wf(dynamic_cast<WaveFunction *>(wf.clone().release())), /*clone returns SamplingFunction ptr*/
        _H(dynamic_cast<Hamiltonian *>(H.clone().release())), /*clone returns ObservableFunction ptr*/
        _nskip_eg(nskip_eg), _blksize_eg(blksize_eg)
{
    if (_wf == nullptr) {
        throw std::runtime_error("[VMC] WaveFunction's clone() did not produce a type derived from WaveFunction."); // check for potential mistakes in implementation
    }
    if (_H == nullptr) {
        throw std::runtime_error("[VMC] Hamiltonian's clone() did not produce a type derived from Hamiltonian."); // check for potential mistakes in implementation
    }
    if (_wf->getTotalNDim() != _H->getTotalNDim()) {
        throw std::invalid_argument("Error VMC: ndim different between wf and H");
    }
    _mci.addSamplingFunction(std::unique_ptr<mci::SamplingFunctionInterface>(_wf));
    _mci.addObservable(std::unique_ptr<mci::ObservableFunctionInterface>(_H), _blksize_eg, _nskip_eg);
}

// --- compute quantities

void VMC::computeEnergy(int Nmc, double * E, double * dE, bool doFindMRT2step, bool doDecorrelation)
{
    MPIVMC::Integrate(_mci, Nmc, E, dE, doFindMRT2step, doDecorrelation);
}
} // namespace vmc
