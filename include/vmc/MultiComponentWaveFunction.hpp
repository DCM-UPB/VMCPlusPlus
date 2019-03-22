#ifndef VMC_MULTICOMPONENTWAVEFUNCTION_HPP
#define VMC_MULTICOMPONENTWAVEFUNCTION_HPP


#include "vmc/WaveFunction.hpp"

#include <vector>



class MultiComponentWaveFunction: public WaveFunction{
private:
    std::vector<WaveFunction *> _wfs;

    mci::SamplingFunctionInterface * _clone() const final {
        auto newwf = new MultiComponentWaveFunction(_nspacedim, _npart, _flag_vd1, _flag_d1vd1, _flag_d2vd1);
        for (auto & wf : _wfs) {
            newwf->addWaveFunction(wf);
        }
        return newwf;
    }
public:
    MultiComponentWaveFunction(const int &nspacedim, const int &npart, bool flag_vd1=false, bool flag_d1vd1=false, bool flag_d2vd1=false):
    WaveFunction(nspacedim, npart, 0, 0, flag_vd1, flag_d1vd1, flag_d2vd1){}
    ~MultiComponentWaveFunction() override{
        _wfs.clear();
    }


    void addWaveFunction(WaveFunction * wf);

    void setVP(const double *vp) override;

    void getVP(double *vp) override;

    void protoFunction(const double * in, double * out) override;

    double acceptanceFunction(const double * protoold, const double * protonew) const override;

    void computeAllDerivatives(const double *x) override;

    double computeWFValue(const double * protovalues) const override;
};


#endif
