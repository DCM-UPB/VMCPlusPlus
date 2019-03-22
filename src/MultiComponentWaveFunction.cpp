#include "vmc/MultiComponentWaveFunction.hpp"

#include <stdexcept>



void MultiComponentWaveFunction::computeAllDerivatives(const double *x){
    for (WaveFunction * wf : _wfs){
        wf->computeAllDerivatives(x);
    }

    // first derivative
    for (int i=0; i<getTotalNDim(); ++i){
        double d1 = 0.;
        for (WaveFunction * wf : _wfs){
            d1 += wf->getD1DivByWF(i);
        }
        _setD1DivByWF(i, d1);
    }
    // second derivative
    for (int i=0; i<getTotalNDim(); ++i){
        double d2 = 0.;
        for (unsigned int iwf=0; iwf<_wfs.size(); ++iwf){
            d2 += _wfs[iwf]->getD2DivByWF(i);
            for (unsigned int jwf=iwf+1; jwf<_wfs.size(); ++jwf){
                d2 += 2. * _wfs[iwf]->getD1DivByWF(i) * _wfs[jwf]->getD1DivByWF(i);
            }
        }
        _setD2DivByWF(i, d2);
    }
    // first variational
    if (hasVD1()){
        int contvp = 0;
        for (WaveFunction * wf : _wfs){
            for (int ivp=0; ivp<wf->getNVP(); ++ivp){
                _setVD1DivByWF(ivp+contvp, wf->getVD1DivByWF(ivp));
            }
            contvp += wf->getNVP();
        }
    }
    // first cross derivative
    if (hasD1VD1()){
        for (int i=0; i<getTotalNDim(); ++i){
            int contvp = 0;
            for (WaveFunction * wf : _wfs){
                for (int ivp=0; ivp<wf->getNVP(); ++ivp){
                    _setD1VD1DivByWF(i, ivp+contvp, wf->getD1VD1DivByWF(i, ivp));
                }
                contvp += wf->getNVP();
            }
            for (unsigned int iwf=0; iwf<_wfs.size(); ++iwf){
                contvp = 0;
                for (unsigned int jwf=0; jwf<_wfs.size(); ++jwf){
                    if (iwf != jwf){
                        for (int ivp=0; ivp<_wfs[jwf]->getNVP(); ++ivp){
                            _setD1VD1DivByWF(i, ivp+contvp, getD1VD1DivByWF(i, ivp+contvp) + _wfs[iwf]->getD1DivByWF(i) * _wfs[jwf]->getVD1DivByWF(ivp));
                        }
                    }
                    contvp += _wfs[jwf]->getNVP();
                }
            }
        }
    }
    // second cross derivative
    if (hasD2VD1()){
        for (int i=0; i<getTotalNDim(); ++i){
             int contvp = 0;
             for ( WaveFunction * wf : _wfs){
                 for (int ivp=0; ivp<wf->getNVP(); ++ivp){
                    _setD2VD1DivByWF(i, ivp+contvp, wf->getD2VD1DivByWF(i, ivp));
                 }
                 contvp += wf->getNVP();
             }

             for (unsigned int iwf=0; iwf<_wfs.size(); ++iwf){
                 contvp = 0;
                 for (unsigned int jwf=0; jwf<_wfs.size(); ++jwf){
                     if (iwf != jwf){
                         for (int ivp=0; ivp<_wfs[jwf]->getNVP(); ++ivp){
                             _setD2VD1DivByWF(i, ivp+contvp, getD2VD1DivByWF(i, ivp+contvp)
                                              + _wfs[iwf]->getD2DivByWF(i) * _wfs[jwf]->getVD1DivByWF(ivp)
                                              + 2. * _wfs[iwf]->getD1DivByWF(i) * _wfs[jwf]->getD1VD1DivByWF(i, ivp)
                                              );
                         }
                     }
                     contvp += _wfs[jwf]->getNVP();
                 }
             }

             for (unsigned int iwf=0; iwf<_wfs.size()-1; ++iwf){
                 for (unsigned int jwf=iwf+1; jwf<_wfs.size(); ++jwf){
                     contvp = 0;
                     for (unsigned int kwf=0; kwf<_wfs.size(); ++kwf){
                        if ((kwf != iwf) && (kwf != jwf)){
                            for (int ivp=0; ivp<_wfs[kwf]->getNVP(); ++ivp){
                                _setD2VD1DivByWF(i, ivp+contvp, getD2VD1DivByWF(i, ivp+contvp) + 2. * _wfs[iwf]->getD1DivByWF(i) * _wfs[jwf]->getD1DivByWF(i) * _wfs[kwf]->getVD1DivByWF(ivp));
                            }
                        }
                        contvp += _wfs[kwf]->getNVP();
                     }
                 }
             }
        }
    }
}

double MultiComponentWaveFunction::computeWFValue(const double * protovalues) const
{
    double out = 1.;
    int contproto = 0;
    for (WaveFunction * wf : _wfs){
        out *= wf->computeWFValue(protovalues+contproto);
        contproto += wf->getNProto();
    }
    return out;
}

double MultiComponentWaveFunction::acceptanceFunction(const double * protoold, const double * protonew) const {
    double acceptance = 1.;
    int contproto = 0;
    for (WaveFunction * wf : _wfs){
        acceptance *= wf->acceptanceFunction(protoold+contproto, protonew+contproto);
        contproto += wf->getNProto();
    }
    return acceptance;
}

void MultiComponentWaveFunction::protoFunction(const double * in, double * out){
    int contproto = 0;
    for (WaveFunction * wf : _wfs){
        wf->protoFunction(in, out+contproto);
        contproto += wf->getNProto();
    }
}


void MultiComponentWaveFunction::getVP(double *vp){
    int contvp = 0;
    for (WaveFunction * wf : _wfs){
        wf->getVP(vp+contvp);
        contvp += wf->getNVP();
    }
}


void MultiComponentWaveFunction::setVP(const double *vp){
    int contvp = 0;
    for (WaveFunction * wf : _wfs){
        wf->setVP(vp+contvp);
        contvp += wf->getNVP();
    }
}


void MultiComponentWaveFunction::addWaveFunction(WaveFunction * wf){
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


void MultiComponentWaveFunction::_newToOld(const mci::WalkerState &wlk) {
    for (auto & wf : _wfs) {
        wf->newToOld(wlk);
    }
}

void MultiComponentWaveFunction::_oldToNew() {
    for (auto & wf : _wfs) {
        wf->oldToNew();
    }
}
