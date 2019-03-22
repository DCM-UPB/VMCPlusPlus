#include "vmc/TwoBodyJastrow.hpp"

#include <cmath>
// #include <iostream>
//
//
// using namespace std;


void TwoBodyJastrow::protoFunction(const double * x, double * protov)
{
    protov[0] = 0.;
    for (int i=0; i<getNPart()-1; ++i){
        for (int j=i+1; j<getNPart(); ++j){
            // cout << "i = " << i << "    j = " << j << endl;
            // cout << "x = " << _pah->getParticleArray(x, i)[0] << "    " << _pah->getParticleArray(x, i)[1] << "    " << _pah->getParticleArray(x, i)[2] << endl;
            // cout << "y = " << _pah->getParticleArray(x, j)[0] << "    " << _pah->getParticleArray(x, j)[1] << "    " << _pah->getParticleArray(x, j)[2] << endl;
            protov[0] += _u2->u(_pah->getParticleArray(x, i), _pah->getParticleArray(x, j));
            // cout << "protov[0] = " << protov[0] << endl;
            // cout << endl;
        }
    }
}

double TwoBodyJastrow::acceptanceFunction(const double * protoold, const double * protonew) const
{
    return exp(2.0 * (protonew[0] - protoold[0]) );   // the factor 2 comes from the fact that the wf must be squared (sampling from psi^2)
}


void TwoBodyJastrow::computeAllDerivatives(const double *x)
{
    double * d1_divbywf = _getD1DivByWF();
    for (int i=0; i<getTotalNDim(); ++i) { d1_divbywf[i] = 0.; }

    double * d2_divbywf = _getD2DivByWF();
    for (int i=0; i<getTotalNDim(); ++i) { d2_divbywf[i] = 0.; }

    double * vd1_divbywf = _getVD1DivByWF();;
    if (hasVD1() || hasD1VD1()){
        for (int i=0; i<getNVP(); ++i) { vd1_divbywf[i] = 0.; }
    }

    double ** d1vd1_divbywf = _getD1VD1DivByWF();
    if (hasD1VD1()){
        for (int i=0; i<getTotalNDim(); ++i) {
            for (int j=0; j<getNVP(); ++j) { d1vd1_divbywf[i][j] = 0.; }
        }
    }

    double ** d2vd1_divbywf = _getD2VD1DivByWF();
    if (hasD2VD1()){
        for (int i=0; i<getTotalNDim(); ++i) {
            for (int j=0; j<getNVP(); ++j) { d2vd1_divbywf[i][j] = 0.; }
        }
    }

    // --- compute the "pure" terms of the derivatives (they will completed with cross terms afterwards)
    for (int i=0; i<getNPart()-1; ++i){
        for (int j=i+1; j<getNPart(); ++j){
            // pre-compute the pseudo-potential derivatives
            _u2->computeAllDerivatives(_pah->getParticleArray(x, i), _pah->getParticleArray(x, j));

            for (int idim=0; idim<getNSpaceDim(); ++idim){
                const int jdim = idim + getNSpaceDim();
                const int ii = idim + i*getNSpaceDim();
                const int ij = idim + j*getNSpaceDim();
                // first derivatives
                d1_divbywf[ii] += _u2->getD1(idim);
                d1_divbywf[ij] += _u2->getD1(jdim);
                // second derivatives
                d2_divbywf[ii] += _u2->getD2(idim);
                d2_divbywf[ij] += _u2->getD2(jdim);
                // first cross derivatives
                if (hasD1VD1()){
                    for (int ivp=0; ivp<getNVP(); ++ivp){
                        d1vd1_divbywf[ii][ivp] += _u2->getD1VD1(idim, ivp);
                        d1vd1_divbywf[ij][ivp] += _u2->getD1VD1(jdim, ivp);
                    }
                }
                // second cross derivatives
                if (hasD2VD1()){
                    for (int ivp=0; ivp<getNVP(); ++ivp){
                        d2vd1_divbywf[ii][ivp] += _u2->getD2VD1(idim, ivp);
                        d2vd1_divbywf[ij][ivp] += _u2->getD2VD1(jdim, ivp);
                    }
                }
            }
            // variational first derivatives
            if (hasVD1()){
                for (int ivp=0; ivp<getNVP(); ++ivp){
                    vd1_divbywf[ivp] += _u2->getVD1(ivp);
                }
            }
        }
    }
    // --- complete the computation of the derivatives
    // second cross derivatives
    if (hasD2VD1()){
        for (int i=0; i<getTotalNDim(); ++i){
            for (int ivp=0; ivp<getNVP(); ++ivp){
                d2vd1_divbywf[i][ivp] += d1_divbywf[i]*d1_divbywf[i]*vd1_divbywf[ivp] + d2_divbywf[i]*vd1_divbywf[ivp] + 2.*d1_divbywf[i]*d1vd1_divbywf[i][ivp];
            }
        }
    }
    // first cross derivative
    if (hasD1VD1()){
        for (int i=0; i<getTotalNDim(); ++i){
            for (int ivp=0; ivp<getNVP(); ++ivp){
                d1vd1_divbywf[i][ivp] += d1_divbywf[i] * vd1_divbywf[ivp];
            }
        }
    }
    // second derivative
    for (int i=0; i<getTotalNDim(); ++i){
        d2_divbywf[i] += d1_divbywf[i]*d1_divbywf[i];
    }
}

double TwoBodyJastrow::computeWFValue(const double * protovalues) const
{
    return exp(2.0 * protovalues[0]);
}
