#ifndef MULTI_COMPONENT_WAVE_FUNCTION
#define MULTI_COMPONENT_WAVE_FUNCTION


#include "WaveFunction.hpp"


#include <vector>
#include <stdexcept>



class MultiComponentWaveFunction: public WaveFunction{
private:
    std::vector<WaveFunction *> _wfs;

public:
    MultiComponentWaveFunction(const int &nspacedim, const int &npart, bool flag_vd1=true, bool flag_d1vd1=true, bool flag_d2vd1=true):
    WaveFunction(nspacedim, npart, 0, 0, flag_vd1, flag_d1vd1, flag_d2vd1){}

    void addWaveFunction(WaveFunction * wf){
        if (wf->getNSpaceDim() != getNSpaceDim()){
            throw std::invalid_argument( "Provided wf's getNSpaceDim() is not valid" );
        }
        if (wf->getNPart() != getNPart()){
            throw std::invalid_argument( "Provided wf's getNPart() is not valid" );
        }
        if (wf->hasVD1() != hasVD1()){
            throw std::invalid_argument( "Provided wf's hasVD1() is not valid" );
        }
        if (wf->hasD1VD1() != hasD1VD1()){
            throw std::invalid_argument( "Provided wf's hasD1VD1() is not valid" );
        }
        if (wf->hasD2VD1() != hasD2VD1()){
            throw std::invalid_argument( "Provided wf's hasD2VD1() is not valid" );
        }

        _wfs.push_back(wf);
        setNProto( getNProto() + wf->getNProto() );
        setNVP( getNVP() + wf->getNVP() );
    }

    void setVP(const double *vp){
        int contvp = 0;
        for (WaveFunction * wf : _wfs){
            wf->setVP(vp+contvp);
            contvp += wf->getNVP();
        }
    }

    void getVP(double *vp){
        int contvp = 0;
        for (WaveFunction * wf : _wfs){
            wf->getVP(vp+contvp);
            contvp += wf->getNVP();
        }
    }

    void samplingFunction(const double * in, double * out){

    }

    double getAcceptance(){
        
    }

    void computeAllDerivatives(const double *x){
        for (WaveFunction * wf : _wfs){
            wf->computeAllDerivatives();
        }
    }

};


#endif
