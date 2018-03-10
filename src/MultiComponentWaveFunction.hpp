#ifndef MULTI_COMPONENT_WAVE_FUNCTION
#define MULTI_COMPONENT_WAVE_FUNCTION


#include "WaveFunction.hpp"

#include <vector>



class MultiComponentWaveFunction: public WaveFunction{
private:
    std::vector<WaveFunction *> _wfs;

public:
    MultiComponentWaveFunction(const int &nspacedim, const int &npart, bool flag_vd1=true, bool flag_d1vd1=true, bool flag_d2vd1=true):
    WaveFunction(nspacedim, npart, 0, 0, flag_vd1, flag_d1vd1, flag_d2vd1){}
    virtual ~MultiComponentWaveFunction(){
        _wfs.clear();
    }


    void addWaveFunction(WaveFunction * wf);

    void setVP(const double *vp);

    void getVP(double *vp);

    void samplingFunction(const double * in, double * out);

    double getAcceptance(const double * protoold, const double * protonew);

    void computeAllDerivatives(const double *x);

};


#endif
