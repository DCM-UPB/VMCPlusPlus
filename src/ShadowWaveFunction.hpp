#ifndef SHADOW_WAVE_FUNCTION
#define SHADOW_WAVE_FUNCTION

#include "WaveFunction.hpp"
#include "PureShadowWaveFunction.hpp"

#include <vector>
#include <random>

/*

    exp( - (R-S)^2 / tau)

*/

class ShadowWaveFunction: public WaveFunction{
private:
    double _tau;
    int _num_swf_sampling;
    double ** _s1;   // sampled shadow1 coordinate
    double ** _s2;   // sampled shadow2 coordinate
    std::vector<PureShadowWaveFunction *> _pswfs;
    std::random_device _rdev;
    std::mt19937_64 _rgen;

public:
    ShadowWaveFunction(const double &tau, const int &num_swf_sampling, const int &nspacedim, const int &npart, bool flag_vd1=true, bool flag_d1vd1=true, bool flag_d2vd1=true);
    ~ShadowWaveFunction();


    // --- manage the pure shadow component
    void addPureShadowWaveFunction(PureShadowWaveFunction * pswf);
    PureShadowWaveFunction * getPureShadowWaveFunction(const int &i){return _pswfs[i];}


    // --- manage the variational parameters
    void getVP(double *vp);
    void setVP(const double *vp);


    // --- sampling-related methods
    void samplingFunction(const double * x, double * proto);
    double getAcceptance(const double * protoold, const double * protonew);


    // --- derivatives
    void computeAllDerivatives(const double *x);

};



#endif
