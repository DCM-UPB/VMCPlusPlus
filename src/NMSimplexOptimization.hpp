#ifndef NM_SIMPLEX_OPTIMIZATION
#define NM_SIMPLEX_OPTIMIZATION

#include "WFOptimization.hpp"
#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "MCIntegrator.hpp"

#include <gsl/gsl_multimin.h>

class NMSimplexOptimization: public WFOptimization
{
protected:
    const long _Nmc;
    const double _iota, _kappa, _lambda;
    const double _rstart, _rend;
public:
    NMSimplexOptimization(WaveFunction * wf, Hamiltonian * H, MCI * mci, const long &Nmc, const double &iota, const double &kappa, const double &lambda, const double &rstart = 1.0, const double &rend = 0.01): WFOptimization(wf, H, mci), _Nmc(Nmc), _iota(iota), _kappa(kappa), _lambda(lambda), _rstart(rstart), _rend(rend) {}
    virtual ~NMSimplexOptimization(){}

    long getNmc(){return _Nmc;}
    double getIota(){return _iota;}
    double getKappa(){return _kappa;}
    double getLambda(){return _lambda;}
    double getRStart(){return _rstart;}
    double getREnd(){return _rend;}

    // optimization
    void optimizeWF();
};

#endif
