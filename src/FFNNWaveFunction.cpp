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
    _ffnn->setInput(in);
    _ffnn->FFPropagate();
    out[0] = pow(_ffnn->getOutput(0), 2);
}


double FFNNWaveFunction::getAcceptance(){
    if ((getProtoOld(0) == 0.) && (getProtoNew(0) != 0.)){
       return 1.;
    } else if ((getProtoOld(0) == 0.) && (getProtoNew(0) == 0.)) {
       return 0.;
    }

    return getProtoNew(0)/getProtoOld(0);
}




// --- wf derivatives

double FFNNWaveFunction::d1(const int &i, const double *in){
    _ffnn->setInput(in);
    _ffnn->FFPropagate();
    return _ffnn->getFirstDerivative(0, i)/_ffnn->getOutput(0);
}


double FFNNWaveFunction::d2(const int &i, const double *in){
    _ffnn->setInput(in);
    _ffnn->FFPropagate();
    return _ffnn->getSecondDerivative(0, i)/_ffnn->getOutput(0);
}


double FFNNWaveFunction::vd1(const int &i, const double *in){
    _ffnn->setInput(in);
    _ffnn->FFPropagate();
    return _ffnn->getVariationalFirstDerivative(0, i)/_ffnn->getOutput(0);
}



// -- Constructor and destructor


FFNNWaveFunction::FFNNWaveFunction(const int &nspacedim, const int &npart, FeedForwardNeuralNetwork * ffnn)
    :WaveFunction(nspacedim, npart, 1, ffnn->getNBeta()){
    if (ffnn->getNInput() != nspacedim*npart)
        throw std::invalid_argument( "FFNN number of inputs does not fit the nspacedime and npart" );

    if (ffnn->getNOutput() != 1)
        throw std::invalid_argument( "FFNN number of output does not fit the wave function requirement (only one value)" );

    _ffnn = ffnn;
}



FFNNWaveFunction::~FFNNWaveFunction(){
    _ffnn = 0;
}
