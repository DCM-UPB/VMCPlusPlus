#include "FFNNWaveFunction.hpp"

#include <cmath>



// --- interface for manipulating the variational parameters
void FFNNWaveFunction::setVP(const double *vp){
    for (int i=0; i<_ffnn->getNBeta(); ++i)
        _ffnn->setBeta(i, vp[i]);
}


void FFNNWaveFunction::getVP(double *vp){
    for (int i=0; i<_ffnn->getNBeta(); ++i)
        vp[i] = _ffnn->getBeta(i);
}




// --- methods herited from MCISamplingFunctionInterface

void FFNNWaveFunction::samplingFunction(const double * in, double * out){
    computeAllInternalValues(in);
    out[0] = pow(getWFValue(), 2);
}


double FFNNWaveFunction::getAcceptance(){
    if ((getProtoOld(0) == 0.) && (getProtoNew(0) != 0.)){
        return 1.;
    } else if ((getProtoOld(0) == 0.) && (getProtoNew(0) == 0.)) {
        return 0.;
    }

    return getProtoNew(0)/getProtoOld(0);
}




// --- computation of the derivatives

void FFNNWaveFunction::computeAllInternalValues(const double *in){
    _ffnn->setInput(in);
    _ffnn->FFPropagate();

    _wf_value = _ffnn->getOutput(0);

    _ffnn->getFirstDerivative(0, _d1_logwf);
    for (int i=0; i<getNDim(); ++i) _d1_logwf[i] /= _wf_value;

    _ffnn->getSecondDerivative(0, _d2_logwf);
    for (int i=0; i<getNDim(); ++i) _d2_logwf[i] /= _wf_value;

    _ffnn->getVariationalFirstDerivative(0, _vd1_logwf);
    for (int i=0; i<getNVP(); ++i) _vd1_logwf[i] /= _wf_value;

}




// -- Constructor and destructor


FFNNWaveFunction::FFNNWaveFunction(const int &nspacedim, const int &npart, FeedForwardNeuralNetwork * ffnn)
    :WaveFunction(nspacedim, npart, 1, ffnn->getNBeta()){
    if (ffnn->getNInput() != nspacedim*npart)
        throw std::invalid_argument( "FFNN number of inputs does not fit the nspacedime and npart" );

    if (ffnn->getNOutput() != 1)
        throw std::invalid_argument( "FFNN number of output does not fit the wave function requirement (only one value)" );

    _d1_logwf = new double[getNDim()];
    _d2_logwf = new double[getNDim()];
    _vd1_logwf = new double[getNVP()];

    _ffnn = ffnn;
}



FFNNWaveFunction::~FFNNWaveFunction(){
    delete[] _d1_logwf;
    delete[] _d2_logwf;
    delete[] _vd1_logwf;

    _ffnn = 0;
}
