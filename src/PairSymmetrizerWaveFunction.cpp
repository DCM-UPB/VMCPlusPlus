#include "PairSymmetrizerWaveFunction.hpp"

#include <cmath>

void PairSymmetrizerWaveFunction::computeAllDerivatives(const double * x)
{
    const int ndim = getTotalNDim();
    double xh[ndim]; // helper array for input
    const double normf = 1. / (_npart*(_npart-1)/2 + 1); // normalizing factor (inverse number of evaluations)
    const double normf2 = -normf; // negative factor for odd permutations in antisym case

    // initialize
    for (int i=0; i<ndim; ++i) xh[i] = x[i];

    // evaluate unswapped wf
    _computeStandardDerivatives(x, normf);

    // add all pair-swapped wfs
    for (int i=0; i<_npart; ++i) {
        for (int j=i+1; j<_npart; ++j) {
            _swapParticles(xh, i, j);

            // evaluate and add swap wf
            if (!_flag_antisymmetric) _addSwapDerivatives(xh, normf);
            else _addSwapDerivatives(xh, normf2);

            _swapParticles(xh, j, i); // swap back
        }
    }
}


void PairSymmetrizerWaveFunction::samplingFunction(const double * in, double * out)
{
    const int ndim = getTotalNDim();
    double outh[_wf->getNProto()], inh[ndim]; // helper arrays for input/output

    const double normf = 1. / (_npart*(_npart-1)/2 + 1); // normalizing factor

    // initialize
    for (int i=0; i<ndim; ++i) inh[i] = in[i];

    // evaluate unswapped wf
    _wf->samplingFunction(in, outh);
    out[0] = normf*_wf->computeWFValue(outh);

    // add all pair-swapped wfs
    for (int i=0; i<_npart; ++i) {
        for (int j=i+1; j<_npart; ++j) {
            _swapParticles(inh, i, j);

            // evaluate and add swap wf
            _wf->samplingFunction(inh, outh);
            if (!_flag_antisymmetric) out[0] += normf*_wf->computeWFValue(outh);
            else out[0] -= normf*_wf->computeWFValue(outh);

            _swapParticles(inh, j, i); // swap back
        }
    }
}
