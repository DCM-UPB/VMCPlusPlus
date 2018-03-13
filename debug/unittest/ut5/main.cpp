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



class PolynomialU2: public TwoBodyPseudoPotential{
/*
    u(r) = a * r^2 + b * r^3
*/

private:
    double _a, _b;

public:
    PolynomialU2(EuclideanMetric * em, double a, double b):
    TwoBodyPseudoPotential(em, 2, true, true, true){
        _a = a;
        _b = b;
    }
    ~PolynomialU2(){}

    void setVP(const double *vp){
        _a=vp[0]; _b=vp[1];
    }
    void getVP(double *vp){
        vp[0]=_a; vp[1]=_b;
    }

    double ur(const double &r){
        return _a * pow(r, 2) + _b * pow(r, 3);
    }

    double urD1(const double &r){
        return 2. * _a * r + 3. * _b * pow(r, 2);
    }

    double urD2(const double &r){
        return 2. * _a + 6. * _b * r;
    }

    void urVD1(const double &r, double * vd1){
        vd1[0] = r*r;
        vd1[1] = r*r*r;
    }

    void urD1VD1(const double &r, double * d1vd1){
        d1vd1[0] = 2.*r;
        d1vd1[1] = 3.*r*r;
    }

    void urD2VD1(const double &r, double * d2vd1){
        d2vd1[0] = 2.;
        d2vd1[1] = 6.*r;
    }
};




int main(){
    using namespace std;

    // constants
    const int NSPACEDIM = 3;
    const int NPART = 2;
    const double DX = 0.0001;
    const double TINY = 0.1;

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd;
    rgen = mt19937_64(rdev());
    rgen.seed(18984687);
    rd = uniform_real_distribution<double>(-0.05, 0.05);

    // distance metric and two body-pseudopotential
    EuclideanMetric * em = new EuclideanMetric(NSPACEDIM);
    PolynomialU2 * u2 = new PolynomialU2(em, -0.3, -0.1);
    TwoBodyJastrow * J = new TwoBodyJastrow(NPART, u2);

    // particles position
    double * x = new double[NPART*NSPACEDIM];

    // pick x from a grid
    const double K = 0.5;
    x[0] = 0.0; x[1] = 0.0; x[2] = K;
    x[3] = K; x[4] = 0.0; x[5] = 0.0;
    // x[6] = 0.0; x[7] = K; x[8] = 0.0;
    // x[9] = 0.0; x[10] = 0.0; x[11] = 0.0;

    // add randomness to x
    for (int i=0; i<NPART; ++i){
        for (int j=0; j<NSPACEDIM; ++j){
            x[i*NSPACEDIM+j] += rd(rgen);
        }
    }
    for (int i=0; i<NSPACEDIM*NPART; ++i) cout << x[i] << "    ";
    cout << endl;

    // variational parameters
    double * vp = new double[J->getNVP()];
    J->getVP(vp);


    // compute all the derivatives analytically
    J->computeAllDerivatives(x);


    // initial wave function
    double f, fdx, fmdx, fdvp, fdxdvp, fmdxdvp;
    J->samplingFunction(x, &f); f = exp(f);
    cout << "f = " << f << endl;


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
        x[i] = origx + DX;
        J->samplingFunction(x, &fdx); fdx = exp(fdx);
        x[i] = origx - DX;
        J->samplingFunction(x, &fmdx); fmdx = exp(fmdx);
        const double numderiv = (fdx - 2.*f + fmdx) / (DX*DX*f);

        // cout << "fdx = " << fdx << "    f = " << f << "    fmdx = " << fmdx << endl;
        // cout << "fdx - 2.*f + fmdx = " << fdx - 2.*f + fmdx << "    DX*DX*f = " << DX*DX*f << endl;
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

            cout << "getD1VD1DivByWF(" << i << ", " << j << ") = " << J->getD1VD1DivByWF(i, j) << endl;
            cout << "numderiv = " << numderiv << endl << endl;
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

            cout << "getD2VD1DivByWF(" << i << ", " << j << ") = " << J->getD2VD1DivByWF(i, j) << endl;
            cout << "numderiv = " << numderiv << endl << endl;
            // assert( abs( (J->getD2VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );

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
