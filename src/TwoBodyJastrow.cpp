#include "TwoBodyJastrow.hpp"

#include <math.h>


void TwoBodyJastrow::samplingFunction(const double * x, double * protov){
    protov[0] = 0.;
    for (int i=0; i<getNPart()-1; ++i){
        for (int j=i+1; j<getNPart(); ++j){
            protov[0] += _u2->u(_pah->getParticleArray(x, i), _pah->getParticleArray(x, j));
        }
    }
}


double TwoBodyJastrow::getAcceptance(){
    return exp(2.0 * (getProtoNew(0) - getProtoOld(0)) );   // the factor 2 comes from the fact that the wf must be squared (sampling from psi^2)
}


void TwoBodyJastrow::computeAllDerivatives(const double *x){

    double * d1_divbywf = _getD1DivByWF();
    for (int i=0; i<getTotalNDim(); ++i) d1_divbywf[i] = 0.;

    double * d2_divbywf = _getD2DivByWF();
    for (int i=0; i<getTotalNDim(); ++i) d2_divbywf[i] = 0.;

    double * vd1_divbywf = _getVD1DivByWF();;
    if (hasVD1() || hasD1VD1()){
        for (int i=0; i<getNVP(); ++i) vd1_divbywf[i] = 0.;
    }

    double ** d1vd1_divbywf = _getD1VD1DivByWF();
    if (hasD1VD1()){
        for (int i=0; i<getTotalNDim(); ++i) for (int j=0; j<getNVP(); ++j) d1vd1_divbywf[i][j] = 0.;
    }

    double ** d2vd1_divbywf = _getD2VD1DivByWF();
    if (hasD2VD1()){
        for (int i=0; i<getTotalNDim(); ++i) for (int j=0; j<getNVP(); ++j) d2vd1_divbywf[i][j] = 0.;
    }

    for (int i=0; i<getNPart()-1; ++i){
        for (int j=i+1; j<getNPart(); ++j){
            // pre-compute the pseudo-potential derivatives
            _u2->computeAllDerivatives(_pah->getParticleArray(x, i), _pah->getParticleArray(x, j));

            // first derivative
            for (int idim=0; idim<getNSpaceDim(); ++idim){
                const int jdim = idim + getNSpaceDim();
                const int ii = idim + i*getNSpaceDim();
                const int ij = idim + j*getNSpaceDim();
                d1_divbywf[ii] += _u2->getD1(idim);
                d1_divbywf[ij] += _u2->getD1(jdim);
            }

            // second derivative
            for (int idim=0; idim<getNSpaceDim(); ++idim){
                const int jdim = idim + getNSpaceDim();
                const int ii = idim + i*getNSpaceDim();
                const int ij = idim + j*getNSpaceDim();
                d2_divbywf[ii] += _u2->getD1(idim) * _u2->getD1(idim) + _u2->getD2(idim);
                d2_divbywf[ij] += _u2->getD1(jdim) * _u2->getD1(jdim) + _u2->getD2(jdim);
            }

            // first variational derivative
            if (hasVD1()){
                for (int ivp=0; ivp<getNVP(); ++ivp) vd1_divbywf[ivp] += _u2->getVD1(ivp);
            }

            // first cross derivative
            if (hasD1VD1()){
                for (int idim=0; idim<getNSpaceDim(); ++idim){
                    const int jdim = idim + getNSpaceDim();
                    const int ii = idim + i*getNSpaceDim();
                    const int ij = idim + j*getNSpaceDim();
                    for (int ivp=0; ivp<getNVP(); ++ivp){
                        d1vd1_divbywf[ii][ivp] += _u2->getVD1(ivp) * _u2->getD1(idim) + _u2->getD1VD1(idim, ivp);
                        d1vd1_divbywf[ij][ivp] += _u2->getVD1(ivp) * _u2->getD1(jdim) + _u2->getD1VD1(jdim, ivp);
                    }
                }
            }

            // second cross derivative
            if (hasD2VD1()){
                for (int idim=0; idim<getNSpaceDim(); ++idim){
                    const int jdim = idim + getNSpaceDim();
                    const int ii = idim + i*getNSpaceDim();
                    const int ij = idim + j*getNSpaceDim();
                    for (int ivp=0; ivp<getNVP(); ++ivp){
                        d2vd1_divbywf[ii][ivp] += _u2->getVD1(ivp) * _u2->getD1(idim) * _u2->getD1(idim) + 2. * _u2->getD1VD1(idim, ivp) * _u2->getD1(idim) + _u2->getVD1(ivp) * _u2->getD2(idim) + _u2->getD2VD1(idim, ivp);
                        d2vd1_divbywf[ij][ivp] += _u2->getVD1(ivp) * _u2->getD1(jdim) * _u2->getD1(jdim) + 2. * _u2->getD1VD1(jdim, ivp) * _u2->getD1(jdim) + _u2->getVD1(ivp) * _u2->getD2(jdim) + _u2->getD2VD1(jdim, ivp);
                    }
                }
            }
        }
    }

}
