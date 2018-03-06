#include "EuclidianParticleDistance.hpp"
#include "TwoBodyPseudoPotential.hpp"

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

};




int main(){
    using namespace std;

    // constants
    const int NSPACEDIM = 3;
    const double DX = 0.001;
    const double TINY = 0.01;
    const int NTEST = 1;

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd1, rd2;
    rgen = std::mt19937_64(rdev());
    rd1 = uniform_real_distribution<double>(-1., -0.1);
    rd2 = uniform_real_distribution<double>(0.1, 1.);

    // distance metric and two body-pseudopotential
    EuclidianParticleDistance * dist = new EuclidianParticleDistance(NSPACEDIM);
    He3u2 * u2 = new He3u2(dist);

    // position of two particles
    double * x = new double[NSPACEDIM];
    double * y = new double[NSPACEDIM];

    // analytical derivative
    double * analderivxy = new double[2*NSPACEDIM];
    double * analderivyx = new double[2*NSPACEDIM];


    // do NTEST tests with random x and y
    for (int k=0; k<NTEST; ++k){
        // pick random x and y
        for (int i=0; i<NSPACEDIM; ++i){
            x[i] = rd1(rgen);
            y[i] = rd2(rgen);
        }



        // --- check the first derivatives

        // compute analytical derivative
        u2->d1(x, y, analderivxy);
        u2->d1(y, x, analderivyx);

        // compute value of the u2
        const double f = u2->u(x, y);

        // derivative in respect to x
        for (int i=0; i<NSPACEDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            const double fdx = u2->u(x, y);
            const double numderiv = (fdx - f)/DX;

            // cout << "analderivxy[" << i << "] = " << analderivxy[i] << endl;
            // cout << "analderivyx[" << i+NSPACEDIM << "] = " << analderivyx[i+NSPACEDIM] << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs(analderivxy[i]-numderiv)/numderiv < TINY );
            assert( abs(analderivyx[i+NSPACEDIM]-numderiv)/numderiv < TINY );

            x[i] = origx;
        }

        // derivative in respect to y
        for (int i=0; i<NSPACEDIM; ++i){
            const double origy = y[i];
            y[i] += DX;
            const double fdy = u2->u(x, y);
            const double numderiv = (fdy - f)/DX;

            // cout << "analderivyx[" << i << "] = " << analderivyx[i] << endl;
            // cout << "analderivxy[" << i+NSPACEDIM << "] = " << analderivxy[i+NSPACEDIM] << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs(analderivyx[i]-numderiv)/numderiv < TINY );
            assert( abs(analderivxy[i+NSPACEDIM]-numderiv)/numderiv < TINY );

            y[i] = origy;
        }



        // -- check the second derivatives

        // compute analytical derivative
        u2->d2(x, y, analderivxy);
        u2->d2(y, x, analderivyx);

        // derivative in respect to x
        for (int i=0; i<NSPACEDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            const double fdx = u2->u(x, y);
            x[i] -= 2.*DX;
            const double fmdx = u2->u(x, y);
            const double numderiv = (fdx - 2.*f + fmdx)/(DX*DX);

            // cout << "analderivxy[" << i << "] = " << analderivxy[i] << endl;
            // cout << "analderivyx[" << i+NSPACEDIM << "] = " << analderivyx[i+NSPACEDIM] << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs(analderivxy[i]-numderiv) < TINY );
            assert( abs(analderivyx[i+NSPACEDIM]-numderiv) < TINY );

            x[i] = origx;
        }

        // derivative in respect to y
        for (int i=0; i<NSPACEDIM; ++i){
            const double origy = y[i];
            y[i] += DX;
            const double fdy = u2->u(x, y);
            y[i] -= 2.*DX;
            const double fmdy = u2->u(x, y);
            const double numderiv = (fdy - 2.*f + fmdy)/(DX*DX);

            // cout << "analderivyx[" << i << "] = " << analderivyx[i] << endl;
            // cout << "analderivxy[" << i+NSPACEDIM << "] = " << analderivxy[i+NSPACEDIM] << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs(analderivyx[i]-numderiv) < TINY );
            assert( abs(analderivxy[i+NSPACEDIM]-numderiv) < TINY );

            y[i] = origy;
        }


    }





    delete[] analderivxy;
    delete[] analderivyx;
    delete[] x;
    delete[] y;
    delete u2;
    delete dist;

    return 0;
}
