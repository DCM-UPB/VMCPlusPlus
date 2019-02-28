#ifndef VMC_MULTICOMPONENTWAVEFUNCTION_HPP
#define VMC_MULTICOMPONENTWAVEFUNCTION_HPP


#include "vmc/WaveFunction.hpp"

#include <vector>



class MultiComponentWaveFunction: public WaveFunction{
private:
    std::vector<WaveFunction *> _wfs;

public:
    MultiComponentWaveFunction(const int &nspacedim, const int &npart, bool flag_vd1=false, bool flag_d1vd1=false, bool flag_d2vd1=false):
    WaveFunction(nspacedim, npart, 0, 0, flag_vd1, flag_d1vd1, flag_d2vd1){}
    ~MultiComponentWaveFunction() override{
        _wfs.clear();
    }


    void addWaveFunction(WaveFunction * wf);

    void setVP(const double *vp) override;

    void getVP(double *vp) override;

    void samplingFunction(const double * in, double * out) override;

    double getAcceptance(const double * protoold, const double * protonew) override;

    void computeAllDerivatives(const double *x) override;

    double computeWFValue(const double * protovalues) override;
};


#endif
