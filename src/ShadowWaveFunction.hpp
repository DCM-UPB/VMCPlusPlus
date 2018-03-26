#ifndef SHADOW_WAVE_FUNCTION
#define SHADOW_WAVE_FUNCTION

/*

    exp( - (R-S)^2 / 4 tau)

*/

class ShadowWaveFunction: public WaveFunction{
public:
    ShadowWaveFunction(const double &tau, const int &nspacedim, const int &npart, bool flag_vd1=true, bool flag_d1vd1=true, bool flag_d2vd1=true):
    WaveFunction(nspacedim, npart, 1, 1, flag_vd1, flag_d1vd1, flag_d2vd1){
        setVP(0, tau);
    }

    virtual ~ShadowWaveFunction(){}

};



#endif
