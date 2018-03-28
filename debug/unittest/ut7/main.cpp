#include "ShadowWaveFunction.hpp"

#include <iostream>
#include <cmath>
#include <assert.h>
#include <stdexcept>


int main(){
    using namespace std;

    const double TAU = 0.5;
    const int NSWFSAMPLING = 1000;

    ShadowWaveFunction * swf;


    // --- check constructor
    bool flag_exception = false;
    try {swf = new ShadowWaveFunction(0., NSWFSAMPLING, 3, 3);}
    catch (exception& e) {flag_exception = true;}
    assert(flag_exception);

    flag_exception = false;
    try {swf = new ShadowWaveFunction(-1., NSWFSAMPLING, 3, 3);}
    catch (exception& e) {flag_exception = true;}
    assert(flag_exception);

    flag_exception = false;
    try {swf = new ShadowWaveFunction(TAU, 0, 3, 3);}
    catch (exception& e) {flag_exception = true;}
    assert(flag_exception);

    flag_exception = false;
    try {swf = new ShadowWaveFunction(TAU, -10, 3, 3);}
    catch (exception& e) {flag_exception = true;}
    assert(flag_exception);

    flag_exception = false;
    try {swf = new ShadowWaveFunction(TAU, NSWFSAMPLING, 3, 3, false, true, false);}
    catch (exception& e) {flag_exception = true;}
    assert(flag_exception);

    flag_exception = false;
    try {swf = new ShadowWaveFunction(TAU, NSWFSAMPLING, 3, 3, false, false, true);}
    catch (exception& e) {flag_exception = true;}
    assert(flag_exception);



    // build the actual SWF we will test
    swf = new ShadowWaveFunction(TAU, NSWFSAMPLING, 3, 3);



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






    delete[] vp;
    delete swf;

    return 0;
}
