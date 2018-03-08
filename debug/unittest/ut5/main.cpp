#include "EuclidianParticleDistance.hpp"
#include "TwoBodyPseudoPotential.hpp"
#include "TwoBodyJastrow.hpp"

#include <assert.h>
#include <math.h>
#include <iostream>
#include <random>



class He3u2: public TwoBodyPseudoPotential{
/*
    u(r) = b/r^5
*/

private:
    double _b;

public:
    He3u2(ParticleDistance * dist):
    TwoBodyPseudoPotential(dist){
        _b = 1;
    }

    int getNVP(){return 1;}
    void setVP(const double *vp){_b=vp[0];}
    void getVP(double *vp){vp[0]=_b;}

    double ur(const double &dist){
        return _b/pow(dist, 5);
    }

    double urD1(const double &dist){
        return -5.*_b/pow(dist, 6);
    }

    double urD2(const double &dist){
        return 30.*_b/pow(dist, 7);
    }

    void urVD1(const double &dist, double * vd1){
        vd1[0] = 1./pow(dist, 5);
    }

};




int main(){
    using namespace std;

    // constants
    const int NSPACEDIM = 3;
    const int NPART = 3;
    const double DX = 0.0001;
    const double TINY = 0.1;

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd;
    rgen = std::mt19937_64(rdev());
    rd = uniform_real_distribution<double>(-0.3, 0.3);

    // distance metric and two body-pseudopotential
    EuclidianParticleDistance * dist = new EuclidianParticleDistance(NSPACEDIM);
    He3u2 * u2 = new He3u2(dist);
    TwoBodyJastrow * J = new TwoBodyJastrow(NPART, u2, true, true, true);

    // particles position
    double * x = new double[NPART*NSPACEDIM];

    // pick x from a grid
    const double K = 2.;
    x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
    x[3] = -K; x[4] = 0.0; x[5] = 0.0;
    x[6] = K; x[7] = 0.0; x[8] = 0.0;
    // x[9] = 0.0; x[10] = -K; x[11] = 0.0;
    // x[12] = 0.0; x[13] = K; x[14] = 0.0;
    // x[15] = 0.0; x[16] = 0.0; x[17] = -K;
    // x[18] = 0.0; x[19] = 0.0; x[20] = K;

    // add randomness to x
    for (int i=0; i<NPART; ++i){
        for (int j=0; j<NSPACEDIM; ++j){
            x[i*NSPACEDIM+j] += rd(rgen);
        }
    }

    // variational parameters
    double * vp = new double[J->getNVP()];
    J->getVP(vp);


    // compute all the derivatives analytically
    J->computeAllDerivatives(x);


    // initial wave function
    double f, fdx, fmdx, fdvp;
    J->samplingFunction(x, &f); f = exp(f);


    // --- check the first derivatives
    for (int i=0; i<NPART*NSPACEDIM; ++i){
        const double origx = x[i];
        x[i] += DX;
        J->samplingFunction(x, &fdx); fdx = exp(fdx);
        const double numderiv = (fdx-f)/(DX*f);

        // cout << "getD1DivByWF(" << i <<") = " << J->getD1DivByWF(i) << endl;
        // cout << "numderiv = " << numderiv << endl << endl;
        assert( abs( (J->getD1DivByWF(i) - numderiv)/numderiv) < TINY );

        x[i] = origx;
    }


    // --- check the second derivatives
    for (int i=0; i<NPART*NSPACEDIM; ++i){
        const double origx = x[i];
        x[i] += DX;
        J->samplingFunction(x, &fdx); fdx = exp(fdx);
        x[i] -= 2.*DX;
        J->samplingFunction(x, &fmdx); fmdx = exp(fmdx);
        const double numderiv = (fdx - 2.*f + fmdx) / (DX*DX*f);

        // cout << "getD2DivByWF(" << i << ") = " << J->getD2DivByWF(i) << endl;
        // cout << "numderiv = " << numderiv << endl << endl;
        assert( abs( (J->getD2DivByWF(i) - numderiv)/numderiv) < TINY );

        x[i] = origx;
    }


    // -- check the first variational derivative
    for (int i=0; i<J->getNVP(); ++i){
        const double origvp = vp[i];
        vp[i] += DX;
        J->setVP(vp);
        J->samplingFunction(x, &fdvp); fdvp = exp(fdvp);
        const double numderiv = (fdvp - f)/(DX*f);

        // cout << "getVD1DivByWF(" << i << ") = " << J->getVD1DivByWF(i) << endl;
        // cout << "numderiv = " << numderiv << endl << endl;
        assert( abs( (J->getVD1DivByWF(i) - numderiv)/numderiv ) < TINY );

        vp[i] = origvp;
        J->setVP(vp);
    }




    delete[] vp;
    delete[] x;
    delete J;
    delete u2;
    delete dist;

    return 0;
}
