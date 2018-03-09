#include "EuclideanMetric.hpp"
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
    He3u2(EuclideanMetric * em):
    TwoBodyPseudoPotential(em, 1, true, true, true){
        _b = -1.;
    }

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

    void urD1VD1(const double &dist, double * d1vd1){
        d1vd1[0] = -5./pow(dist, 6);
    }

    void urD2VD1(const double &dist, double * d1vd1){
        d1vd1[0] = 30./pow(dist, 7);
    }

};




int main(){
    using namespace std;

    // constants
    const int NSPACEDIM = 3;
    const int NPART = 3;
    const double DX = 0.001;
    const double TINY = 0.1;

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd;
    rgen = std::mt19937_64(rdev());
    rd = uniform_real_distribution<double>(-0.3, 0.3);

    // distance metric and two body-pseudopotential
    EuclideanMetric * em = new EuclideanMetric(NSPACEDIM);
    He3u2 * u2 = new He3u2(em);
    TwoBodyJastrow * J = new TwoBodyJastrow(NPART, u2, true, true, true);

    // particles position
    double * x = new double[NPART*NSPACEDIM];

    // pick x from a grid
    const double K = 2.;
    x[0] = 0.0; x[1] = 0.0; x[2] = K;
    x[3] = K; x[4] = 0.0; x[5] = 0.0;
    x[6] = 0.0; x[7] = K; x[8] = 0.0;

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
    double f, fdx, fmdx, fdvp, fdxdvp, fmdxdvp;
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


    // --- check the first cross derivative
    for (int i=0; i<NPART*NSPACEDIM; ++i){
        for (int j=0; j<J->getNVP(); ++j){
            const double origx = x[i];
            const double origvp = vp[j];

            x[i] += DX;
            J->samplingFunction(x, &fdx); fdx = exp(fdx);

            x[i] = origx;
            vp[j] += DX;
            J->setVP(vp);
            J->samplingFunction(x, &fdvp); fdvp = exp(fdvp);

            x[i] += DX;
            J->samplingFunction(x, &fdxdvp); fdxdvp = exp(fdxdvp);

            const double numderiv = (fdxdvp - fdx - fdvp + f)/(DX*DX*f);

            // cout << "getD1VD1DivByWF(" << i << ", " << j << ") = " << J->getD1VD1DivByWF(i, j) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs( (J->getD1VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );

            x[i] = origx;
            vp[j] = origvp;
            J->setVP(vp);
        }
    }


    // --- check the second cross derivative
    for (int i=0; i<NPART*NSPACEDIM; ++i){
        for (int j=0; j<J->getNVP(); ++j){
            const double origx = x[i];
            const double origvp = vp[j];

            x[i] = origx + DX;
            J->samplingFunction(x, &fdx); fdx = exp(fdx);

            vp[j] = origvp + DX;
            J->setVP(vp);
            J->samplingFunction(x, &fdxdvp); fdxdvp = exp(fdxdvp);

            x[i] = origx;
            J->samplingFunction(x, &fdvp); fdvp = exp(fdvp);

            x[i] = origx - DX;
            J->samplingFunction(x, &fmdxdvp); fmdxdvp = exp(fmdxdvp);

            vp[j] = origvp;
            J->setVP(vp);
            J->samplingFunction(x, &fmdx); fmdx = exp(fmdx);

            const double numderiv = (fdxdvp - 2.*fdvp + fmdxdvp - fdx + 2.*f - fmdx)/(DX*DX*DX*f);

            // cout << "getD2VD1DivByWF(" << i << ", " << j << ") = " << J->getD2VD1DivByWF(i, j) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs( (J->getD2VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );

            x[i] = origx;
            vp[j] = origvp;
            J->setVP(vp);
        }
    }




    delete[] vp;
    delete[] x;
    delete J;
    delete u2;
    delete em;

    return 0;
}
