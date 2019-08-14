#ifndef VMC_MULTICOMPONENTWAVEFUNCTION_HPP
#define VMC_MULTICOMPONENTWAVEFUNCTION_HPP


#include "vmc/WaveFunction.hpp"

#include <vector>

namespace vmc
{

class MultiComponentWaveFunction final: public WaveFunction
{
private:
    std::vector<WaveFunction *> _wfs;

    mci::SamplingFunctionInterface * _clone() const final
    {
        auto newwf = new MultiComponentWaveFunction(_nspacedim, _npart, _flag_vd1, _flag_d1vd1, _flag_d2vd1);
        for (auto &wf : _wfs) {
            newwf->addWaveFunction(wf);
        }
        return newwf;
    }

    // we contain ProtoFunctionInterfaces as members, so we need to implement these:
    void _newToOld() final;
    void _oldToNew() final;

public:
    MultiComponentWaveFunction(int nspacedim, int npart, bool flag_vd1 = false, bool flag_d1vd1 = false, bool flag_d2vd1 = false):
            WaveFunction(nspacedim, npart, 0, 0, flag_vd1, flag_d1vd1, flag_d2vd1) {}
    ~MultiComponentWaveFunction() final
    {
        _wfs.clear();
    }


    void addWaveFunction(WaveFunction * wf);

    void setVP(const double vp[]) final;

    void getVP(double vp[]) const final;

    void protoFunction(const double in[], double out[]) final;

    double acceptanceFunction(const double protoold[], const double protonew[]) const final;

    double updatedAcceptance(const mci::WalkerState &wlk, const double protoold[], double protonew[]) final;

    void computeAllDerivatives(const double x[]) final;

    double computeWFValue(const double protovalues[]) const final;
};
} // namespace vmc

#endif
