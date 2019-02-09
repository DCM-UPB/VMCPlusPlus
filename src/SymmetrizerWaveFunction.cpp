#include "vmc/SymmetrizerWaveFunction.hpp"

#include <stdexcept>
#include <cmath>

unsigned long SymmetrizerWaveFunction::_npart_factorial()
{
    unsigned long fac = 1;
    const unsigned long npart = static_cast<unsigned long>(_npart);
    for (unsigned long i=2; i<=npart; ++i) fac *= i;
    return fac;
}

void SymmetrizerWaveFunction::_swapPositions(double * x, const int &i, const int &j)
{
    double xh;
    // particle swap (of positions)
    for (int k=0; k<_nspacedim; ++k) {
        xh = x[i*_nspacedim+k];
        x[i*_nspacedim+k] = x[j*_nspacedim+k];

        x[j*_nspacedim+k] = xh;
    }
}

void SymmetrizerWaveFunction::_swapIndices(int * ids, const int &i, const int &j)
{
    int idh;
    // particle swap (of indices)
    for (int k=0; k<_nspacedim; ++k) {
        idh = ids[i*_nspacedim+k];
        ids[i*_nspacedim+k] = ids[j*_nspacedim+k];

        ids[j*_nspacedim+k] = idh;
    }
}

void SymmetrizerWaveFunction::_computeStandardDerivatives(const double * x, const double &normf)
{
    double protov[_wf->getNProto()];
    _wf->samplingFunction(x, protov);
    const double wfvalue = _wf->computeWFValue(protov);
    const double normf2 = normf * wfvalue;

    _wf->computeAllDerivatives(x);

    const int ndim = getTotalNDim();
    for (int ip=0; ip<ndim; ++ip) {
        _setD1DivByWF(ip, normf2*_wf->getD1DivByWF(ip));
        _setD2DivByWF(ip, normf2*_wf->getD2DivByWF(ip));
    }

    if (hasVD1()) {
        for (int ivp=0; ivp<_nvp; ++ivp) {
            _setVD1DivByWF(ivp, normf2*_wf->getVD1DivByWF(ivp));
        }
    }

    if (hasD1VD1()){
        for (int ip=0; ip<ndim; ++ip) {
            for (int ivp=0; ivp<_nvp; ++ivp) {
                _setD1VD1DivByWF(ip, ivp, normf2*_wf->getD1VD1DivByWF(ip, ivp));
            }
        }
    }

    if (hasD2VD1()){
        for (int ip=0; ip<ndim; ++ip){
            for (int ivp=0; ivp<_nvp; ++ivp) {
                _setD2VD1DivByWF(ip, ivp, normf2*_wf->getD2VD1DivByWF(ip, ivp));
            }
        }
    }
}

void SymmetrizerWaveFunction::_addSwapDerivatives(const double * x, const double &normf, const int * ids)
{
    double protov[_wf->getNProto()];
    _wf->samplingFunction(x, protov);
    const double wfvalue = _wf->computeWFValue(protov);
    const double normf2 = normf * wfvalue;

    _wf->computeAllDerivatives(x);

    const int ndim = getTotalNDim();
    for (int ip=0; ip<ndim; ++ip) {
        _setD1DivByWF(ids[ip], getD1DivByWF(ids[ip]) + normf2*_wf->getD1DivByWF(ip));
        _setD2DivByWF(ids[ip], getD2DivByWF(ids[ip]) + normf2*_wf->getD2DivByWF(ip));
    }

    if (hasVD1()) {
        for (int ivp=0; ivp<_nvp; ++ivp) {
            _setVD1DivByWF(ivp, getVD1DivByWF(ivp) + normf2*_wf->getVD1DivByWF(ivp));
        }
    }

    if (hasD1VD1()){
        for (int ip=0; ip<ndim; ++ip) {
            for (int ivp=0; ivp<_nvp; ++ivp) {
                _setD1VD1DivByWF(ids[ip], ivp, getD1VD1DivByWF(ids[ip], ivp) + normf2*_wf->getD1VD1DivByWF(ip, ivp));
            }
        }
    }

    if (hasD2VD1()){
        for (int ip=0; ip<ndim; ++ip){
            for (int ivp=0; ivp<_nvp; ++ivp) {
                _setD2VD1DivByWF(ids[ip], ivp, getD2VD1DivByWF(ids[ip], ivp) + normf2*_wf->getD2VD1DivByWF(ip, ivp));
            }
        }
    }
}



void SymmetrizerWaveFunction::computeAllDerivatives(const double * x)
{
    const int ndim = getTotalNDim();
    double xh[ndim]; // helper array for positions
    int idh[ndim]; // helper array for indices
    int counts[_npart], iter; // counters for heaps algorithm
    bool isOdd = false; // flip for permutation parity

    double protov[_wf->getNProto()];
    samplingFunction(x, protov);
    const double normf = 1. / (_npart_factorial()*computeWFValue(protov)); // normalizing factor
    const double normf2 = -normf; // negative factor for odd permutations in antisym case

    // initialize
    iter = 0;
    for (int i=0; i<_npart; ++i) counts[i] = 0;
    for (int i=0; i<ndim; ++i) {
        xh[i] = x[i];
        idh[i] = i;
    }

    // evaluate unswapped wf
    _computeStandardDerivatives(x, normf);

    // add swapped wfs by heaps algorithm
    while (iter < _npart) {
        if (counts[iter] < iter) {
            if (iter % 2 == 0) {
                _swapPositions(xh, 0, iter);
                _swapIndices(idh, 0, iter);
                isOdd = false;
            }
            else {
                _swapPositions(xh, counts[iter], iter);
                _swapIndices(idh, counts[iter], iter);
                isOdd = true;
            }

            // evaluate and add swap wf
            if (!_flag_antisymmetric || !isOdd) _addSwapDerivatives(xh, normf, idh);
            else _addSwapDerivatives(xh, normf2, idh);

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
                _swapPositions(inh, 0, iter);
                isOdd = false;
            }
            else {
                _swapPositions(inh, counts[iter], iter);
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
