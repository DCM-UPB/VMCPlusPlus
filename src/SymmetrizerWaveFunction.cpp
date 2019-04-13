#include "vmc/SymmetrizerWaveFunction.hpp"

#include <cmath>

namespace vmc
{

unsigned long SymmetrizerWaveFunction::_npart_factorial()
{
    unsigned long fac = 1;
    const auto npart = static_cast<unsigned long>(_npart);
    for (unsigned long i = 2; i <= npart; ++i) { fac *= i; }

    return fac;
}

void SymmetrizerWaveFunction::_swapPositions(double * x, const int i, const int j)
{
    // particle swap (of positions)
    std::swap_ranges(x + i*_nspacedim, x + (i + 1)*_nspacedim, x + j*_nspacedim);
}

void SymmetrizerWaveFunction::_swapIndices(int * ids, const int i, const int j)
{
    // particle swap (of indices)
    std::swap_ranges(ids + i*_nspacedim, ids + (i + 1)*_nspacedim, ids + j*_nspacedim);
}

void SymmetrizerWaveFunction::_computeStandardDerivatives(const double * x, const double normf)
{
    double protov[_wf->getNProto()];
    _wf->protoFunction(x, protov);
    const double wfvalue = _wf->computeWFValue(protov);
    const double normf2 = normf*wfvalue;

    _wf->computeAllDerivatives(x);

    const int ndim = getTotalNDim();
    for (int ip = 0; ip < ndim; ++ip) {
        _setD1DivByWF(ip, normf2*_wf->getD1DivByWF(ip));
        _setD2DivByWF(ip, normf2*_wf->getD2DivByWF(ip));
    }

    if (hasVD1()) {
        for (int ivp = 0; ivp < _nvp; ++ivp) {
            _setVD1DivByWF(ivp, normf2*_wf->getVD1DivByWF(ivp));
        }
    }

    if (hasD1VD1()) {
        for (int ip = 0; ip < ndim; ++ip) {
            for (int ivp = 0; ivp < _nvp; ++ivp) {
                _setD1VD1DivByWF(ip, ivp, normf2*_wf->getD1VD1DivByWF(ip, ivp));
            }
        }
    }

    if (hasD2VD1()) {
        for (int ip = 0; ip < ndim; ++ip) {
            for (int ivp = 0; ivp < _nvp; ++ivp) {
                _setD2VD1DivByWF(ip, ivp, normf2*_wf->getD2VD1DivByWF(ip, ivp));
            }
        }
    }
}

void SymmetrizerWaveFunction::_addSwapDerivatives(const double * x, const double normf, const int * ids)
{
    double protov[_wf->getNProto()];
    _wf->protoFunction(x, protov);
    const double wfvalue = _wf->computeWFValue(protov);
    const double normf2 = normf*wfvalue;

    _wf->computeAllDerivatives(x);

    const int ndim = getTotalNDim();
    for (int ip = 0; ip < ndim; ++ip) {
        _setD1DivByWF(ids[ip], getD1DivByWF(ids[ip]) + normf2*_wf->getD1DivByWF(ip));
        _setD2DivByWF(ids[ip], getD2DivByWF(ids[ip]) + normf2*_wf->getD2DivByWF(ip));
    }

    if (hasVD1()) {
        for (int ivp = 0; ivp < _nvp; ++ivp) {
            _setVD1DivByWF(ivp, getVD1DivByWF(ivp) + normf2*_wf->getVD1DivByWF(ivp));
        }
    }

    if (hasD1VD1()) {
        for (int ip = 0; ip < ndim; ++ip) {
            for (int ivp = 0; ivp < _nvp; ++ivp) {
                _setD1VD1DivByWF(ids[ip], ivp, getD1VD1DivByWF(ids[ip], ivp) + normf2*_wf->getD1VD1DivByWF(ip, ivp));
            }
        }
    }

    if (hasD2VD1()) {
        for (int ip = 0; ip < ndim; ++ip) {
            for (int ivp = 0; ivp < _nvp; ++ivp) {
                _setD2VD1DivByWF(ids[ip], ivp, getD2VD1DivByWF(ids[ip], ivp) + normf2*_wf->getD2VD1DivByWF(ip, ivp));
            }
        }
    }
}

void SymmetrizerWaveFunction::_newToOld()
{
    _wf->newToOld();
}

void SymmetrizerWaveFunction::computeAllDerivatives(const double * x)
{
    const int ndim = getTotalNDim();
    double xh[ndim]; // helper array for positions
    int idh[ndim]; // helper array for indices
    int counts[_npart], iter; // counters for heaps algorithm
    bool isOdd = false; // flip for permutation parity

    double protov[_wf->getNProto()];
    protoFunction(x, protov);
    const double normf = 1./(_npart_factorial()*computeWFValue(protov)); // normalizing factor
    const double normf2 = -normf; // negative factor for odd permutations in antisym case

    // initialize
    iter = 0;
    std::fill(counts, counts + _npart, 0.);
    std::copy(x, x + ndim, xh);
    std::iota(idh, idh + ndim, 0); // range 0..ndim-1

    // evaluate unswapped wf
    _computeStandardDerivatives(x, normf);

    // add swapped wfs by heaps algorithm
    while (iter < _npart) {
        if (counts[iter] < iter) {
            if (iter%2 == 0) {
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
            if (!_flag_antisymmetric || !isOdd) {
                _addSwapDerivatives(xh, normf, idh);
            }
            else {
                _addSwapDerivatives(xh, normf2, idh);
            }

            ++counts[iter];
            iter = 0;
        }
        else {
            counts[iter] = 0;
            ++iter;
        }
    }
}

double SymmetrizerWaveFunction::computeWFValue(const double * protovalues) const
{
    return protovalues[0]; // sign is important here so we calculate the unsquared wf in protoFunction
}

double SymmetrizerWaveFunction::acceptanceFunction(const double * protoold, const double * protonew) const
{
    const double po2 = protoold[0]*protoold[0];
    const double pn2 = protonew[0]*protonew[0];
    if (po2 == 0 && pn2 != 0) { return 1.; }
    if (po2 != 0 && pn2 == 0) { return 0.; }
    return pn2/po2;
}


void SymmetrizerWaveFunction::protoFunction(const double * in, double * out)
{
    const int ndim = getTotalNDim();
    double outh[_wf->getNProto()], inh[ndim]; // helper arrays for input/output
    int counts[_npart], iter; // counters for heaps algorithm
    bool isOdd = false; // flip for permutation parity
    const double normf = 1./_npart_factorial(); // normalizing factor

    // initialize
    iter = 0;
    std::fill(counts, counts + _npart, 0.);
    std::copy(in, in + ndim, inh);

    // evaluate unswapped wf
    _wf->protoFunction(in, outh);
    out[0] = normf*_wf->computeWFValue(outh);

    // add swapped wfs by heaps algorithm
    while (iter < _npart) {
        if (counts[iter] < iter) {
            if (iter%2 == 0) {
                _swapPositions(inh, 0, iter);
                isOdd = false;
            }
            else {
                _swapPositions(inh, counts[iter], iter);
                isOdd = true;
            }

            // evaluate and add swap wf
            _wf->protoFunction(inh, outh);
            if (!_flag_antisymmetric || !isOdd) {
                out[0] += normf*_wf->computeWFValue(outh);
            }
            else {
                out[0] -= normf*_wf->computeWFValue(outh);
            }

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
} // namespace vmc