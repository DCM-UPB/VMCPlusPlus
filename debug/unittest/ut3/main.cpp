#include "EuclidianParticleDistance.hpp"

#include <assert.h>
#include <math.h>
#include <iostream>
#include <random>


int main(){
    using namespace std;
    const int NSPACEDIM = 3;
    const double DX = 0.005;
    const double TINY = 0.01;
    const int NTEST = 10;

    EuclidianParticleDistance * epd = new EuclidianParticleDistance(NSPACEDIM);

    double * x = new double[NSPACEDIM];
    double * y = new double[NSPACEDIM];

    for (int i=0; i<NSPACEDIM; ++i){
        x[i] = 2.;
        y[i] = -1.;
    }



    // --- check that the distance is computed correctly
    assert( epd->dist(x, y) == 3.*sqrt(NSPACEDIM) );
    assert( epd->dist(y, x) == 3.*sqrt(NSPACEDIM) );



    // --- check derivatives
    // array that will store the analytical derivatives
    double * analderivxy = new double[2*NSPACEDIM];
    double * analderivyx = new double[2*NSPACEDIM];

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd1, rd2;
    rgen = std::mt19937_64(rdev());
    rd1 = uniform_real_distribution<double>(-5,-1);
    rd2 = uniform_real_distribution<double>(1,5);

    // make NTEST tests, with random x and y
    for (int k=0; k<NTEST; ++k){
        // pick random x and y
        for (int i=0; i<NSPACEDIM; ++i){
            x[i] = rd1(rgen);
            y[i] = rd2(rgen);
        }

        // compute the distance xy
        const double f = epd->dist(x, y);



        // --- check the first derivative

        epd->distD1(x, y, analderivxy);
        epd->distD1(y, x, analderivyx);

        // check derivative in respect to x
        for (int i=0; i<NSPACEDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            const double fdx = epd->dist(x, y);
            const double numderiv = (fdx - f)/DX;

            // cout << "analderivxy[" <<  i << "] = " << analderivxy[i] << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs(analderivxy[i]-numderiv)/numderiv < TINY );
            assert( abs(analderivyx[i+NSPACEDIM]-numderiv)/numderiv < TINY );

            x[i] = origx;
        }

        // check derivative in respect to y
        for (int i=0; i<NSPACEDIM; ++i){
            const double origy = y[i];
            y[i] += DX;
            const double fdy = epd->dist(x, y);
            const double numderiv = (fdy - f)/DX;

            // cout << "analderivyx[" <<  i << "] = " << analderivyx[i] << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs(analderivyx[i]-numderiv)/numderiv < TINY );
            assert( abs(analderivxy[i+NSPACEDIM]-numderiv)/numderiv < TINY );

            y[i] = origy;
        }



        // --- check the second derivative

        epd->distD2(x, y, analderivxy);
        epd->distD2(y, x, analderivyx);

        // check derivative in respect to x
        for (int i=0; i<NSPACEDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            const double fdx = epd->dist(x, y);
            x[i] -= 2.*DX;
            const double fmdx = epd->dist(x, y);
            const double numderiv = (fdx - 2.*f + fmdx)/(DX*DX);

            // cout << "analderivxy[" <<  i << "] = " << analderivxy[i] << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs(analderivxy[i]-numderiv)/numderiv < TINY );
            assert( abs(analderivyx[i+NSPACEDIM]-numderiv)/numderiv < TINY );

            x[i] = origx;
        }

        // check derivative in respect to y
        for (int i=0; i<NSPACEDIM; ++i){
            const double origy = y[i];
            y[i] += DX;
            const double fdy = epd->dist(x, y);
            y[i] -= 2.*DX;
            const double fmdy = epd->dist(x, y);
            const double numderiv = (fdy - 2.*f + fmdy)/(DX*DX);

            // cout << "analderivyx[" <<  i << "] = " << analderivyx[i] << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs(analderivyx[i]-numderiv)/numderiv < TINY );
            assert( abs(analderivxy[i+NSPACEDIM]-numderiv)/numderiv < TINY );

            y[i] = origy;
        }
    }





    // --- free resources
    delete[] analderivxy;
    delete[] analderivyx;
    delete[] y;
    delete[] x;
    delete epd;

    return 0;
}
