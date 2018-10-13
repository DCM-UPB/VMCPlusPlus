#include "SymmetrizerWaveFunction.hpp"

#include <stdexcept>
#include <cmath>

unsigned long SymmetrizerWaveFunction::_npart_factorial()
{
    unsigned long fac = 1;
    const unsigned long npart = static_cast<unsigned long>(_npart);
    for (unsigned long i=2; i<=npart; ++i) fac *= i;
    return fac;
}

void SymmetrizerWaveFunction::_swapParticles(double * x, const int &i, const int &j)
{
    double xh;
    // particle swap
    for (int k=0; k<_nspacedim; ++k) {
        xh = x[i*_nspacedim+k];
        x[i*_nspacedim+k] = x[j*_nspacedim+k];
        x[j*_nspacedim+k] = xh;
    }
}

void SymmetrizerWaveFunction::_computeStandardDerivatives(const double * x, const double &normf)
{
    _wf->computeAllDerivatives(x);

    const int ndim = getTotalNDim();
    for (int i=0; i<ndim; ++i) {
        _setD1DivByWF(i, normf*_wf->getD1DivByWF(i));
        _setD2DivByWF(i, normf*_wf->getD2DivByWF(i));
    }

    if (hasVD1()) {
        for (int ivp=0; ivp<_nvp; ++ivp) {
            _setVD1DivByWF(ivp, normf*_wf->getVD1DivByWF(ivp));
        }
    }

    if (hasD1VD1()){
        for (int i=0; i<ndim; ++i) {
            for (int ivp=0; ivp<_nvp; ++ivp) {
                _setD1VD1DivByWF(i, ivp, normf*_wf->getD1VD1DivByWF(i, ivp));
            }
        }
    }

    if (hasD2VD1()){
        for (int i=0; i<ndim; ++i){
            for (int ivp=0; ivp<_nvp; ++ivp) {
                _setD2VD1DivByWF(i, ivp, normf*_wf->getD2VD1DivByWF(i, ivp));
            }
        }
    }
}

void SymmetrizerWaveFunction::_addSwapDerivatives(const double * x, const double &normf)
{
    _wf->computeAllDerivatives(x);

    const int ndim = getTotalNDim();
    for (int i=0; i<ndim; ++i) {
        _setD1DivByWF(i, getD1DivByWF(i) + normf*_wf->getD1DivByWF(i));
        _setD2DivByWF(i, getD2DivByWF(i) + normf*_wf->getD2DivByWF(i));
    }

    if (hasVD1()) {
        for (int ivp=0; ivp<_nvp; ++ivp) {
            _setVD1DivByWF(ivp, getVD1DivByWF(ivp) + normf*_wf->getVD1DivByWF(ivp));
        }
    }

    if (hasD1VD1()){
        for (int i=0; i<ndim; ++i) {
            for (int ivp=0; ivp<_nvp; ++ivp) {
                _setD1VD1DivByWF(i, ivp, getD1VD1DivByWF(i, ivp) + normf*_wf->getD1VD1DivByWF(i, ivp));
            }
        }
    }

    if (hasD2VD1()){
        for (int i=0; i<ndim; ++i){
            for (int ivp=0; ivp<_nvp; ++ivp) {
                _setD2VD1DivByWF(i, ivp, getD2VD1DivByWF(i, ivp) + normf*_wf->getD2VD1DivByWF(i, ivp));
            }
        }
    }
}



void SymmetrizerWaveFunction::computeAllDerivatives(const double * x)
{
    const int ndim = getTotalNDim();
    double xh[ndim]; // helper arrays for input/output
    int counts[_npart], iter; // counters for heaps algorithm
    bool isOdd = false; // flip for permutation parity
    const double normf = 1. / _npart_factorial(); // normalizing factor
    const double normf2 = -normf; // negative factor for odd permutations in antisym case

    // initialize
    iter = 0;
    for (int i=0; i<_npart; ++i) counts[i] = 0;
    for (int i=0; i<ndim; ++i) xh[i] = x[i];

    // evaluate unswapped wf
    _computeStandardDerivatives(x, normf);

    // add swapped wfs by heaps algorithm
    while (iter < _npart) {
        if (counts[iter] < iter) {
            if (iter % 2 == 0) {
                _swapParticles(xh, 0, iter);
                isOdd = false;
            }
            else {
                _swapParticles(xh, counts[iter], iter);
                isOdd = true;
            }

            // evaluate and add swap wf
            if (!_flag_antisymmetric || !isOdd) _addSwapDerivatives(xh, normf);
            else _addSwapDerivatives(xh, normf2);

            ++counts[iter];
            iter = 0;
        }
        else {
            counts[iter] = 0;
            ++iter;
        }
    }
}

double SymmetrizerWaveFunction::computeWFValue(const double * protovalues)
{
    return protovalues[0]; // sign is important here so we calculate the unsquared wf in samplingFunction
}

double SymmetrizerWaveFunction::getAcceptance(const double * protoold, const double * protonew)
{
    double po2 = protoold[0]*protoold[0], pn2 = protonew[0]*protonew[0];
    if (po2 == 0 && pn2 != 0) return 1.;
    else if (po2 != 0 && pn2 == 0) return 0.;
    else return pn2 / po2;
}


void SymmetrizerWaveFunction::samplingFunction(const double * in, double * out)
{
    const int ndim = getTotalNDim();
    double outh[_wf->getNProto()], inh[ndim]; // helper arrays for input/output
    int counts[_npart], iter; // counters for heaps algorithm
    bool isOdd = false; // flip for permutation parity
    const double normf = 1. / _npart_factorial(); // normalizing factor

    // initialize
    iter = 0;
    for (int i=0; i<_npart; ++i) counts[i] = 0.;
    for (int i=0; i<ndim; ++i) inh[i] = in[i];

    // evaluate unswapped wf
    _wf->samplingFunction(in, outh);
    out[0] = normf*_wf->computeWFValue(outh);

    // add swapped wfs by heaps algorithm
    while (iter < _npart) {
        if (counts[iter] < iter) {
            if (iter % 2 == 0) {
                _swapParticles(inh, 0, iter);
                isOdd = false;
            }
            else {
                _swapParticles(inh, counts[iter], iter);
                isOdd = true;
            }

            // evaluate and add swap wf
            _wf->samplingFunction(inh, outh);
            if (!_flag_antisymmetric || !isOdd) out[0] += normf*_wf->computeWFValue(outh);
            else out[0] -= normf*_wf->computeWFValue(outh);

            ++counts[iter];
            iter = 0;
        }
        else {
            counts[iter] = 0;
            ++iter;
        }
    }
}


void SymmetrizerWaveFunction::getVP(double * vp)
{
    _wf->getVP(vp);
}


void SymmetrizerWaveFunction::setVP(const double * vp)
{
    _wf->setVP(vp);
}
