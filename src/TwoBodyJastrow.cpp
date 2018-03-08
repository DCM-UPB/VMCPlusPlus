#include "TwoBodyJastrow.hpp"

#include <math.h>


void TwoBodyJastrow::samplingFunction(const double * x, double * protov){
    protov[0] = 0.;
    for (int i=0; i<getNPart()-1; ++i){
        for (int j=i+1; j<getNPart(); ++j){
            protov[0] += _u2->u(_pam->getParticleArray(x, i), _pam->getParticleArray(x, j));
        }
    }
}


double TwoBodyJastrow::getAcceptance(){
    return exp(2.0 * (getProtoNew(0) - getProtoOld(0)) );   // the factor 2 comes from the fact that the wf must be squared (sampling from psi^2)
}


void TwoBodyJastrow::computeAllDerivatives(const double *x){

    double * d1_divbywf = getD1DivByWF();
    for (int i=0; i<getTotalNDim(); ++i) d1_divbywf[i] = 0.;

    double * d2_divbywf = getD2DivByWF();
    for (int i=0; i<getTotalNDim(); ++i) d2_divbywf[i] = 0.;



    double * ud1 = new double[2*getNSpaceDim()];
    double * ud1_2 = new double[2*getNSpaceDim()];
    double * ud2 = new double[2*getNSpaceDim()];

    for (int i=0; i<getNPart()-1; ++i){
        for (int j=i+1; j<getNPart(); ++j){

            _u2->d1(_pam->getParticleArray(x, i), _pam->getParticleArray(x, j), ud1);
            for (int i=0; i<2*getNSpaceDim(); ++i) ud1_2[i] = ud1[i]*ud1[i];

            _u2->d2(_pam->getParticleArray(x, i), _pam->getParticleArray(x, j), ud2);


            _pam->addArrayToParticleArray(d1_divbywf, i, ud1);
            _pam->addArrayToParticleArray(d1_divbywf, j, ud1+getNSpaceDim());

            _pam->addArrayToParticleArray(d2_divbywf, i, ud1_2);
            _pam->addArrayToParticleArray(d2_divbywf, j, ud1_2+getNSpaceDim());

            _pam->addArrayToParticleArray(d2_divbywf, i, ud2);
            _pam->addArrayToParticleArray(d2_divbywf, j, ud2+getNSpaceDim());
        }
    }


    delete[] ud2;
    delete[] ud1_2;
    delete[] ud1;
}
