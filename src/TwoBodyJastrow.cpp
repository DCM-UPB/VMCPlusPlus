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

    double * vd1_divbywf = getVD1DivByWF();;
    if (hasVD1() || hasD1VD1()){
        for (int i=0; i<getNVP(); ++i) vd1_divbywf[i] = 0.;
    }

    double ** d1vd1_divbywf = getD1VD1DivByWF();
    if (hasD1VD1()){
        for (int i=0; i<getTotalNDim(); ++i) for (int j=0; j<getNVP(); ++j) d1vd1_divbywf[i][j] = 0.;
    }

    double ** d2vd1_divbywf = getD2VD1DivByWF();
    if (hasD2VD1()){
        for (int i=0; i<getTotalNDim(); ++i) for (int j=0; j<getNVP(); ++j) d2vd1_divbywf[i][j] = 0.;
    }



    double * ud1 = new double[2*getNSpaceDim()];
    double * ud1_2 = new double[2*getNSpaceDim()];
    double * ud2 = new double[2*getNSpaceDim()];
    double * uvd1 = new double[getNVP()];
    double ** ud1vd1 = new double*[2*getNSpaceDim()];
    for (int i=0; i<2*getNSpaceDim(); ++i) ud1vd1[i] = new double[getNVP()];
    double ** ud2vd1 = new double*[2*getNSpaceDim()];
    for (int i=0; i<2*getNSpaceDim(); ++i) ud2vd1[i] = new double[getNVP()];

    for (int i=0; i<getNPart()-1; ++i){
        for (int j=i+1; j<getNPart(); ++j){
            // compute factors used for Jastrow derivatives
            _u2->d1(_pam->getParticleArray(x, i), _pam->getParticleArray(x, j), ud1);
            for (int i=0; i<2*getNSpaceDim(); ++i) ud1_2[i] = ud1[i]*ud1[i];

            _u2->d2(_pam->getParticleArray(x, i), _pam->getParticleArray(x, j), ud2);

            if (hasVD1() || hasD1VD1() || hasD2VD1()){
                _u2->vd1(_pam->getParticleArray(x, i), _pam->getParticleArray(x, j), uvd1);
            }

            if (hasD1VD1() || hasD2VD1()){
                _u2->d1vd1(_pam->getParticleArray(x, i), _pam->getParticleArray(x, j), ud1vd1);
            }

            if (hasD2VD1()){
                _u2->d2vd1(_pam->getParticleArray(x, i), _pam->getParticleArray(x, j), ud2vd1);
            }


            // update derivatives
            _pam->addArrayToParticleArray(d1_divbywf, i, ud1);
            _pam->addArrayToParticleArray(d1_divbywf, j, ud1+getNSpaceDim());

            _pam->addArrayToParticleArray(d2_divbywf, i, ud1_2);
            _pam->addArrayToParticleArray(d2_divbywf, j, ud1_2+getNSpaceDim());

            _pam->addArrayToParticleArray(d2_divbywf, i, ud2);
            _pam->addArrayToParticleArray(d2_divbywf, j, ud2+getNSpaceDim());

            if (hasVD1()){
                for (int ivd1=0; ivd1<getNVP(); ++ivd1) vd1_divbywf[ivd1] += uvd1[ivd1];
            }

            if (hasD1VD1()){
                for (int idim=0; idim<getNSpaceDim(); ++idim){
                    const int ii = idim + i*getNSpaceDim();
                    const int ij = idim + j*getNSpaceDim();
                    for (int ivd1=0; ivd1<getNVP(); ++ivd1){
                        d1vd1_divbywf[ii][ivd1] += uvd1[ivd1] * ud1[idim] + ud1vd1[idim][ivd1];
                        d1vd1_divbywf[ij][ivd1] += uvd1[ivd1] * ud1[idim+getNSpaceDim()] + ud1vd1[idim+getNSpaceDim()][ivd1];
                    }
                }
            }

            if (hasD2VD1()){
                for (int idim=0; idim<getNSpaceDim(); ++idim){
                    const int ii = idim + i*getNSpaceDim();
                    const int ij = idim + j*getNSpaceDim();
                    for (int ivd1=0; ivd1<getNVP(); ++ivd1){
                        d2vd1_divbywf[ii][ivd1] += uvd1[ivd1] * ud1[idim] * ud1[idim] + 2. * ud1vd1[idim][ivd1] * ud1[idim] + uvd1[ivd1] * ud2[idim] + ud2vd1[idim][ivd1];
                        d2vd1_divbywf[ij][ivd1] += uvd1[ivd1] * ud1[idim+getNSpaceDim()] * ud1[idim+getNSpaceDim()] + 2. * ud1vd1[idim+getNSpaceDim()][ivd1] * ud1[idim+getNSpaceDim()] + uvd1[ivd1] * ud2[idim+getNSpaceDim()] + ud2vd1[idim+getNSpaceDim()][ivd1];
                    }
                }
            }
        }
    }


    for (int i=0; i<2*getNSpaceDim(); ++i) delete[] ud2vd1[i];
    delete[] ud2vd1;
    for (int i=0; i<2*getNSpaceDim(); ++i) delete[] ud1vd1[i];
    delete[] ud1vd1;
    delete[] uvd1;
    delete[] ud2;
    delete[] ud1_2;
    delete[] ud1;
}
