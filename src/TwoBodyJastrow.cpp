#include "TwoBodyJastrow.hpp"

#include <math.h>


void TwoBodyJastrow::samplingFunction(const double * x, double * protov){
    protov[0] = 0.;
    for (int i=0; i<getNPart()-1; ++i){
        for (int j=i+1; j<getNPart(); ++j){
            protov[0] -= _u2->u(_ptc_ary->getParticleArray(x, i), _ptc_ary->getParticleArray(x, j));
        }
    }
}


double TwoBodyJastrow::getAcceptance(){
    return exp(getProtoNew(0) - getProtoOld(0));
}


void TwoBodyJastrow::computeAllDerivatives(const double *x){

    double * d1_divbywf = getD1DivByWF();
    for (int i=0; i<getTotalNDim(); ++i) d1_divbywf[i] = 0;

    double * vec = new double[2*getNSpaceDim()];
    for (int i=0; i<getNPart()-1; ++i){
        for (int j=i+1; j<getNPart(); ++j){
            _u2->d1(_ptc_ary->getParticleArray(x, i), _ptc_ary->getParticleArray(x, j), vec);

            for (int i=0; i<2*getNSpaceDim(); ++i){
                vec[i] = -vec[i];
            }

            _ptc_ary->addArrayToParticleArray(d1_divbywf, i, vec);
            _ptc_ary->addArrayToParticleArray(d1_divbywf, j, vec+getNSpaceDim());
        }
    }



    double * d2_divbywf = getD2DivByWF();
    for (int i=0; i<getTotalNDim(); ++i) d2_divbywf[i] = 0;

    for (int i=0; i<getNPart()-1; ++i){
        for (int j=i+1; j<getNPart(); ++j){
            _u2->d2(_ptc_ary->getParticleArray(x, i), _ptc_ary->getParticleArray(x, j), vec);

            for (int i=0; i<2*getNSpaceDim(); ++i){
                vec[i] = -vec[i];
            }

            _ptc_ary->addArrayToParticleArray(d2_divbywf, i, vec);
            _ptc_ary->addArrayToParticleArray(d2_divbywf, j, vec+getNSpaceDim());
        }
    }



    delete[] vec;

}
