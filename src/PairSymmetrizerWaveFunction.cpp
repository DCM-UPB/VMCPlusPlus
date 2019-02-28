#include "vmc/PairSymmetrizerWaveFunction.hpp"

#include <cmath>

void PairSymmetrizerWaveFunction::computeAllDerivatives(const double * x)
{
    const int ndim = getTotalNDim();
    double xh[ndim]; // helper array for positions
    int idh[ndim]; // helper array for indices

    double protov[_wf->getNProto()];
    samplingFunction(x, protov);
    const double normf = 1. / ((_npart*(_npart-1)/2 + 1)*computeWFValue(protov)); // normalizing factor
    const double normf2 = -normf; // negative factor for odd permutations in antisym case

    // initialize
    for (int i=0; i<ndim; ++i) {
        xh[i] = x[i];
        idh[i] = i;
    }

    // evaluate unswapped wf
    _computeStandardDerivatives(x, normf);

    // add all pair-swapped wfs
    for (int i=0; i<_npart; ++i) {
        for (int j=i+1; j<_npart; ++j) {
            _swapPositions(xh, i, j);
            _swapIndices(idh, i, j);

            // evaluate and add swap wf
            if (!_flag_antisymmetric) { _addSwapDerivatives(xh, normf, idh);
            } else { _addSwapDerivatives(xh, normf2, idh); }

            _swapPositions(xh, j, i); // swap back
            _swapIndices(idh, j, i);
        }
    }
}


void PairSymmetrizerWaveFunction::samplingFunction(const double * in, double * out)
{
    const int ndim = getTotalNDim();
    double outh[_wf->getNProto()], inh[ndim]; // helper arrays for input/output

    const double normf = 1. / (_npart*(_npart-1)/2 + 1); // normalizing factor

    // initialize
    for (int i=0; i<ndim; ++i) { inh[i] = in[i]; }

    // evaluate unswapped wf
    _wf->samplingFunction(in, outh);
    out[0] = normf*_wf->computeWFValue(outh);

    // add all pair-swapped wfs
    for (int i=0; i<_npart; ++i) {
        for (int j=i+1; j<_npart; ++j) {
            _swapPositions(inh, i, j);

            // evaluate and add swap wf
            _wf->samplingFunction(inh, outh);
            if (!_flag_antisymmetric) { out[0] += normf*_wf->computeWFValue(outh);
            } else { out[0] -= normf*_wf->computeWFValue(outh); }

            _swapPositions(inh, j, i); // swap back
        }
    }
}
