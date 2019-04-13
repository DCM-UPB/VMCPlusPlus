#ifndef VMC_NMSIMPLEXOPTIMIZATION_HPP
#define VMC_NMSIMPLEXOPTIMIZATION_HPP

#include "mci/MCIntegrator.hpp"
#include "vmc/Hamiltonian.hpp"
#include "vmc/WFOptimization.hpp"
#include "vmc/WaveFunction.hpp"

namespace vmc
{

class NMSimplexOptimization: public WFOptimization
{
protected:
    const int _Nmc;
    const double _iota, _kappa, _lambda;
    const double _rstart, _rend;
    const size_t _max_n_iter;
public:
    NMSimplexOptimization(WaveFunction * wf, Hamiltonian * H, mci::MCI * mci, int Nmc, double iota, double kappa, double lambda, double rstart = 1.0, double rend = 0.01, size_t max_n_iter = 0):
            WFOptimization(wf, H, mci), _Nmc(Nmc), _iota(iota), _kappa(kappa), _lambda(lambda), _rstart(rstart), _rend(rend), _max_n_iter(static_cast<size_t>(max_n_iter)) {}
    ~NMSimplexOptimization() override = default;

    int getNmc() { return _Nmc; }
    double getIota() { return _iota; }
    double getKappa() { return _kappa; }
    double getLambda() { return _lambda; }
    double getRStart() { return _rstart; }
    double getREnd() { return _rend; }
    size_t getMaxNIter() { return _max_n_iter; }

    // optimization
    void optimizeWF() override;
};
} // namespace vmc

#endif
