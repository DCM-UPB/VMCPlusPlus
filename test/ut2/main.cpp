#include "vmc/EuclideanMetric.hpp"

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

    EuclideanMetric em(NSPACEDIM);

    double x[NSPACEDIM];
    double y[NSPACEDIM];

    for (int i=0; i<NSPACEDIM; ++i){
        x[i] = 2.;
        y[i] = -1.;
    }



    // --- check that the distance is computed correctly
    assert( em.dist(x, y) == 3.*sqrt(NSPACEDIM) );
    assert( em.dist(y, x) == 3.*sqrt(NSPACEDIM) );



    // --- check derivatives
    // array that will store the analytical derivatives
    double analderivxy[2*NSPACEDIM];
    double analderivyx[2*NSPACEDIM];

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd1, rd2;
    rgen = std::mt19937_64(rdev());
    rgen.seed(18984687);
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
        const double f = em.dist(x, y);



        // --- check the first derivative

        em.distD1(x, y, analderivxy);
        em.distD1(y, x, analderivyx);

        // check derivative in respect to x
        for (int i=0; i<NSPACEDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            const double fdx = em.dist(x, y);
            const double numderiv = (fdx - f)/DX;

            // cout << "analderivxy[" <<  i << "] = " << analderivxy[i] << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs(analderivxy[i]-numderiv)/numderiv < TINY );
            assert( fabs(analderivyx[i+NSPACEDIM]-numderiv)/numderiv < TINY );

            x[i] = origx;
        }

        // check derivative in respect to y
        for (int i=0; i<NSPACEDIM; ++i){
            const double origy = y[i];
            y[i] += DX;
            const double fdy = em.dist(x, y);
            const double numderiv = (fdy - f)/DX;

            // cout << "analderivyx[" <<  i << "] = " << analderivyx[i] << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs(analderivyx[i]-numderiv)/numderiv < TINY );
            assert( fabs(analderivxy[i+NSPACEDIM]-numderiv)/numderiv < TINY );

            y[i] = origy;
        }



        // --- check the second derivative

        em.distD2(x, y, analderivxy);
        em.distD2(y, x, analderivyx);

        // check derivative in respect to x
        for (int i=0; i<NSPACEDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            const double fdx = em.dist(x, y);
            x[i] -= 2.*DX;
            const double fmdx = em.dist(x, y);
            const double numderiv = (fdx - 2.*f + fmdx)/(DX*DX);

            // cout << "analderivxy[" <<  i << "] = " << analderivxy[i] << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs(analderivxy[i]-numderiv)/numderiv < TINY );
            assert( fabs(analderivyx[i+NSPACEDIM]-numderiv)/numderiv < TINY );

            x[i] = origx;
        }

        // check derivative in respect to y
        for (int i=0; i<NSPACEDIM; ++i){
            const double origy = y[i];
            y[i] += DX;
            const double fdy = em.dist(x, y);
            y[i] -= 2.*DX;
            const double fmdy = em.dist(x, y);
            const double numderiv = (fdy - 2.*f + fmdy)/(DX*DX);

            // cout << "analderivyx[" <<  i << "] = " << analderivyx[i] << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs(analderivyx[i]-numderiv)/numderiv < TINY );
            assert( fabs(analderivxy[i+NSPACEDIM]-numderiv)/numderiv < TINY );

            y[i] = origy;
        }
    }

    return 0;
}
