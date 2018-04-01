#include "ShadowWaveFunction.hpp"

#include <iostream>
#include <cmath>
#include <assert.h>
#include <stdexcept>


int main(){
    using namespace std;

    const double TAU = 0.5;
    const int NSWFSAMPLING = 1000000;
    const int NDIM = 3;
    const int NPART = 3;
    const double TINY = 0.1;
    const long SEED = 12507153;

    ShadowWaveFunction * swf;


    // --- check constructor
    bool flag_exception = false;
    try {swf = new ShadowWaveFunction(0., NSWFSAMPLING, NDIM, NPART);}
    catch (exception& e) {flag_exception = true;}
    assert(flag_exception);

    flag_exception = false;
    try {swf = new ShadowWaveFunction(-1., NSWFSAMPLING, NDIM, NPART);}
    catch (exception& e) {flag_exception = true;}
    assert(flag_exception);

    flag_exception = false;
    try {swf = new ShadowWaveFunction(TAU, 0, NDIM, NPART);}
    catch (exception& e) {flag_exception = true;}
    assert(flag_exception);

    flag_exception = false;
    try {swf = new ShadowWaveFunction(TAU, -10, NDIM, NPART);}
    catch (exception& e) {flag_exception = true;}
    assert(flag_exception);

    flag_exception = false;
    try {swf = new ShadowWaveFunction(TAU, NSWFSAMPLING, NDIM, NPART, false, true, false);}
    catch (exception& e) {flag_exception = true;}
    assert(flag_exception);

    flag_exception = false;
    try {swf = new ShadowWaveFunction(TAU, NSWFSAMPLING, NDIM, NPART, false, false, true);}
    catch (exception& e) {flag_exception = true;}
    assert(flag_exception);



    // build the actual SWF we will test
    swf = new ShadowWaveFunction(TAU, NSWFSAMPLING, NDIM, NPART);
    swf->setSeed(SEED);



    // --- check variational parameters behaviour
    assert( swf->getNVP() == 1 );
    double * vp = new double[1];
    swf->getVP(vp);
    assert( vp[0] == TAU );
    vp[0] = 7.;
    swf->setVP(vp);
    swf->getVP(vp);
    assert( vp[0] == 7. );
    vp[0] = TAU;
    swf->setVP(vp);



    // --- check sampling and acceptance
    double * x = new double[swf->getTotalNDim()];
    x[0] = 1.; x[1] = 0.; x[2] = 0.;
    x[3] = 0.; x[4] = 1.; x[5] = 0.;
    x[6] = 0.; x[7] = 0.; x[8] = 1.;
    double * protoold = new double[1];
    swf->samplingFunction(x, protoold);
    assert( protoold[0] == 1. );
    x[0] = 0.; x[1] = 3.; x[2] = -1.25;
    x[3] = 1.; x[4] = 0.; x[5] = 0.5;
    x[6] = 2.; x[7] = 1.; x[8] = -1.5;
    double * protonew = new double[1];
    swf->samplingFunction(x, protonew);
    assert( protonew[0] == 1. );
    assert( swf->getAcceptance(protoold, protonew) == 1. );



    // --- check derivatives
    swf->computeAllDerivatives(x);

    // d1
    // cout << " --- Check D1:" << endl;
    for (int i=0; i<swf->getTotalNDim(); ++i){
        // cout << "swf->getD1DivByWF(" << i << ") = " << swf->getD1DivByWF(i) << endl;
        assert( abs(swf->getD1DivByWF(i)) < TINY);
    }

    // d2
    // cout << " --- Check D2:" << endl;
    for (int i=0; i<swf->getTotalNDim(); ++i){
        // cout << "swf->getD2DivByWF(" << i << ") = " << swf->getD2DivByWF(i) << endl;
        assert( abs(swf->getD2DivByWF(i)) < TINY);
    }

    // vd1
    cout << " --- Check VD1:" << endl;
    for (int i=0; i<swf->getNVP(); ++i){
        cout << "swf->getVD1DivByWF(" << i << ") = " << swf->getVD1DivByWF(i) << endl;
        cout << "correct vd1 = " << NDIM*NPART/(2.*TAU) << endl;
        assert( abs(swf->getVD1DivByWF(i) - NDIM*NPART/(2.*TAU)) < TINY);
    }

    // d1vd1
    // cout << " --- Check D1VD1:" << endl;
    for (int i=0; i<swf->getTotalNDim(); ++i){
        // cout << "swf->getD1VD1DivByWF(" << i << ", 0) = " << swf->getD1VD1DivByWF(i, 0) << endl;
        assert( abs(swf->getD1VD1DivByWF(i, 0)) < TINY);
    }

    // d2vd1
    cout << " --- Check D2VD1:" << endl;
    for (int i=0; i<swf->getTotalNDim(); ++i){
        cout << "swf->getD2VD1DivByWF(" << i << ", 0) = " << swf->getD2VD1DivByWF(i, 0) << endl;
        // assert( abs(swf->getD2VD1DivByWF(i, 0)) < TINY);
    }





    delete[] protoold;
    delete[] protonew;
    delete[] x;
    delete[] vp;
    delete swf;

    return 0;
}
