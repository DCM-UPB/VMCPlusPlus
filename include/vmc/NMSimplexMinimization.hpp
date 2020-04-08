#ifndef VMC_NMSIMPLEXMINIMIZATION_HPP
#define VMC_NMSIMPLEXMINIMIZATION_HPP

#include "vmc/VMC.hpp"

namespace vmc
{

class NMSimplexMinimization
{
protected:
    const int64_t _Nmc;
    const double _iota, _kappa, _lambda;
    const double _rstart, _rend;
    const size_t _max_n_iter;
public:
    NMSimplexMinimization(int64_t Nmc, double iota, double kappa, double lambda, double rstart = 1.0, double rend = 0.01, size_t max_n_iter = 0):
            _Nmc(Nmc), _iota(iota), _kappa(kappa), _lambda(lambda), _rstart(rstart), _rend(rend), _max_n_iter(static_cast<size_t>(max_n_iter)) {}

    int64_t getNmc() { return _Nmc; }
    double getIota() { return _iota; }
    double getKappa() { return _kappa; }
    double getLambda() { return _lambda; }
    double getRStart() { return _rstart; }
    double getREnd() { return _rend; }
    size_t getMaxNIter() { return _max_n_iter; }

    // optimization
    void minimizeEnergy(VMC &vmc);
};
} // namespace vmc

#endif
