#include "TwoBodyJastrow.hpp"



void TwoBodyJastrow::samplingFunction(const double * x, double * protov){
    protov[0] = 0.;
    for (int i=0; i<getNPart()-1; ++i){
        for (int j=i+1; j<getNPart(); ++j){
            protov[0] -= _u2->u(_ppos->getParticlePosition(x, i), _ppos->getParticlePosition(x, j));
        }
    }
}


double TwoBodyJastrow::getAcceptance(){
    return exp(getProtoNew(0) - getProtoOld(0));
}


void TwoBodyJastrow::computeAllDerivatives(const double *x){

    double * d1_

}
