#include "vmc/EuclideanMetric.hpp"
#include "vmc/TwoBodyPseudoPotential.hpp"

#include <assert.h>
#include <math.h>
#include <iostream>
#include <random>

#include "TestVMCFunctions.hpp"


int main(){
    using namespace std;

    // constants
    const int NSPACEDIM = 3;
    const double DX = 0.001;
    const double TINY = 0.01;
    const int NTEST = 10;

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd;
    rgen = mt19937_64(rdev());
    rgen.seed(18984687);
    rd = uniform_real_distribution<double>(0.25, -0.25);

    // distance metric and two body-pseudopotential
    EuclideanMetric * em = new EuclideanMetric(NSPACEDIM);
    He3u2 * u2 = new He3u2(em);
    PolynomialU2 * poly_u2 = new PolynomialU2(em, -0.3, -0.1);

    // position of two particles
    double x[NSPACEDIM];
    double y[NSPACEDIM];

    // variational parameters
    double vp[u2->getNVP()];
    u2->getVP(vp);
    double poly_vp[poly_u2->getNVP()];
    poly_u2->getVP(poly_vp);

    // pick x from a grid
    const double K = 0.75;
    x[0] = 0.0; x[1] = 0.0; x[2] = K;
    y[0] = K; y[1] = 0.0; y[2] = 0.0;

    // do NTEST tests with random x and y
    for (int k=0; k<NTEST; ++k){
        // pick random x and y
        for (int i=0; i<NSPACEDIM; ++i){
            x[i] += rd(rgen);
            y[i] += rd(rgen);
        }


        // compute value of the u2
        const double f = u2->u(x, y);
        const double pf = poly_u2->u(x, y);


        // compute all the analytical derivative
        u2->computeAllDerivatives(x, y);
        poly_u2->computeAllDerivatives(x, y);



        // --- check the first derivatives

        // derivative in respect to x
        for (int i=0; i<NSPACEDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            const double fdx = u2->u(x, y);
            const double numderiv = (fdx - f)/DX;

            // cout << "u2->getD1(" << i << ") = " << u2->getD1(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs((u2->getD1(i)-numderiv)/numderiv) < TINY );

            x[i] = origx;
        }
        for (int i=0; i<NSPACEDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            const double fdx = poly_u2->u(x, y);
            const double numderiv = (fdx - pf)/DX;

            // cout << "poly_u2->getD1(" << i << ") = " << poly_u2->getD1(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs((poly_u2->getD1(i)-numderiv)/numderiv) < TINY );

            x[i] = origx;
        }

        // derivative in respect to y
        for (int i=0; i<NSPACEDIM; ++i){
            const double origy = y[i];
            y[i] += DX;
            const double fdy = u2->u(x, y);
            const double numderiv = (fdy - f)/DX;

            // cout << "u2->getD1(" << i+NSPACEDIM << ") = " << u2->getD1(i+NSPACEDIM) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs((u2->getD1(i+NSPACEDIM)-numderiv)/numderiv) < TINY );

            y[i] = origy;
        }
        for (int i=0; i<NSPACEDIM; ++i){
            const double origy = y[i];
            y[i] += DX;
            const double fdy = poly_u2->u(x, y);
            const double numderiv = (fdy - pf)/DX;

            // cout << "poly_u2->getD1(" << i+NSPACEDIM << ") = " << poly_u2->getD1(i+NSPACEDIM) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs((poly_u2->getD1(i+NSPACEDIM)-numderiv)/numderiv) < TINY );

            y[i] = origy;
        }



        // -- check the second derivatives

        // derivative in respect to x
        for (int i=0; i<NSPACEDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            const double fdx = u2->u(x, y);
            x[i] -= 2.*DX;
            const double fmdx = u2->u(x, y);
            const double numderiv = (fdx - 2.*f + fmdx)/(DX*DX);

            // cout << "u2->getD2(" << i << ") = " << u2->getD2(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs( (u2->getD2(i)-numderiv)/numderiv ) < TINY );

            x[i] = origx;
        }
        for (int i=0; i<NSPACEDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            const double fdx = poly_u2->u(x, y);
            x[i] -= 2.*DX;
            const double fmdx = poly_u2->u(x, y);
            const double numderiv = (fdx - 2.*pf + fmdx)/(DX*DX);

            // cout << "poly_u2->getD2(" << i << ") = " << poly_u2->getD2(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs( (poly_u2->getD2(i)-numderiv)/numderiv ) < TINY );

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

            // cout << "u2->getD2(" << i+NSPACEDIM << ") = " << u2->getD2(i+NSPACEDIM) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs( (u2->getD2(i+NSPACEDIM)-numderiv)/numderiv ) < TINY );

            y[i] = origy;
        }
        for (int i=0; i<NSPACEDIM; ++i){
            const double origy = y[i];
            y[i] += DX;
            const double fdy = poly_u2->u(x, y);
            y[i] -= 2.*DX;
            const double fmdy = poly_u2->u(x, y);
            const double numderiv = (fdy - 2.*pf + fmdy)/(DX*DX);

            // cout << "poly_u2->getD2(" << i+NSPACEDIM << ") = " << poly_u2->getD2(i+NSPACEDIM) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs( (poly_u2->getD2(i+NSPACEDIM)-numderiv)/numderiv ) < TINY );

            y[i] = origy;
        }



        // -- check the first variational derivative

        for (int i=0; i<u2->getNVP(); ++i){
            const double origvp = vp[i];
            vp[i] = origvp + DX;
            u2->setVP(vp);
            const double fdp = u2->u(x,y);
            const double numderiv = (fdp - f)/DX;

            // cout << "u2->getVD1(" << i << ") = " << u2->getVD1(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs(( u2->getVD1(i)-numderiv)/numderiv ) < TINY );

            vp[i] = origvp;
            u2->setVP(vp);
        }
        for (int i=0; i<poly_u2->getNVP(); ++i){
            const double origvp = poly_vp[i];
            poly_vp[i] = origvp + DX;
            poly_u2->setVP(poly_vp);
            const double fdp = poly_u2->u(x,y);
            const double numderiv = (fdp - pf)/DX;

            // cout << "poly_u2->getVD1(" << i << ") = " << poly_u2->getVD1(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs(( poly_u2->getVD1(i)-numderiv)/numderiv ) < TINY );

            poly_vp[i] = origvp;
            poly_u2->setVP(poly_vp);
        }



        // --- check the first cross derivative

        // derivative in respect to x
        for (int i=0; i<NSPACEDIM; ++i){
            for (int j=0; j<u2->getNVP(); ++j){
                const double origx = x[i];
                const double origvp = vp[j];

                x[i] += DX;
                const double fdx = u2->u(x, y);

                x[i] = origx;
                vp[j] += DX;
                u2->setVP(vp);
                const double fdp = u2->u(x, y);

                x[i] += DX;
                const double fdxdp = u2->u(x, y);

                const double numderiv = (fdxdp - fdx - fdp + f)/(DX*DX);

                // cout << "u2->getD1VD1(" << i << ", " << j << ") = " << u2->getD1VD1(i, j) << endl;
                // cout << "numderiv = " << numderiv << endl << endl;
                assert( fabs( (u2->getD1VD1(i, j)-numderiv)/numderiv ) < TINY );

                x[i] = origx;
                vp[j] = origvp;
                u2->setVP(vp);
            }
        }
        for (int i=0; i<NSPACEDIM; ++i){
            for (int j=0; j<poly_u2->getNVP(); ++j){
                const double origx = x[i];
                const double origvp = poly_vp[j];

                x[i] += DX;
                const double fdx = poly_u2->u(x, y);

                x[i] = origx;
                poly_vp[j] += DX;
                poly_u2->setVP(poly_vp);
                const double fdp = poly_u2->u(x, y);

                x[i] += DX;
                const double fdxdp = poly_u2->u(x, y);

                const double numderiv = (fdxdp - fdx - fdp + pf)/(DX*DX);

                // cout << "poly_u2->getD1VD1(" << i << ", " << j << ") = " << poly_u2->getD1VD1(i, j) << endl;
                // cout << "numderiv = " << numderiv << endl << endl;
                assert( fabs( (poly_u2->getD1VD1(i, j)-numderiv)/numderiv ) < TINY );

                x[i] = origx;
                poly_vp[j] = origvp;
                poly_u2->setVP(poly_vp);
            }
        }

        // derivative in respect to y
        for (int i=0; i<NSPACEDIM; ++i){
            for (int j=0; j<u2->getNVP(); ++j){
                const double origy = y[i];
                const double origvp = vp[j];

                y[i] += DX;
                const double fdx = u2->u(x, y);

                y[i] = origy;
                vp[j] += DX;
                u2->setVP(vp);
                const double fdp = u2->u(x, y);

                y[i] += DX;
                const double fdxdp = u2->u(x, y);

                const double numderiv = (fdxdp - fdx - fdp + f)/(DX*DX);

                // cout << "u2->getD1VD1(" << i+NSPACEDIM << ", " << j << ") = " << u2->getD1VD1(i+NSPACEDIM, j) << endl;
                // cout << "numderiv = " << numderiv << endl << endl;
                assert( fabs( (u2->getD1VD1(i+NSPACEDIM, j)-numderiv)/numderiv ) < TINY );

                y[i] = origy;
                vp[j] = origvp;
                u2->setVP(vp);
            }
        }
        for (int i=0; i<NSPACEDIM; ++i){
            for (int j=0; j<poly_u2->getNVP(); ++j){
                const double origy = y[i];
                const double origvp = poly_vp[j];

                y[i] += DX;
                const double fdx = poly_u2->u(x, y);

                y[i] = origy;
                poly_vp[j] += DX;
                poly_u2->setVP(poly_vp);
                const double fdp = poly_u2->u(x, y);

                y[i] += DX;
                const double fdxdp = poly_u2->u(x, y);

                const double numderiv = (fdxdp - fdx - fdp + pf)/(DX*DX);

                // cout << "poly_u2->getD1VD1(" << i+NSPACEDIM << ", " << j << ") = " << poly_u2->getD1VD1(i+NSPACEDIM, j) << endl;
                // cout << "numderiv = " << numderiv << endl << endl;
                assert( fabs( (poly_u2->getD1VD1(i+NSPACEDIM, j)-numderiv)/numderiv ) < TINY );

                y[i] = origy;
                poly_vp[j] = origvp;
                poly_u2->setVP(poly_vp);
            }
        }



        // --- check the second cross derivative

        // derivative in respect to x
        for (int i=0; i<NSPACEDIM; ++i){
            for (int j=0; j<u2->getNVP(); ++j){
                const double origx = x[i];
                const double origvp = vp[j];

                x[i] = origx + DX;
                const double fdx = u2->u(x, y);

                x[i] = origx;
                vp[j] += DX;
                u2->setVP(vp);
                const double fdp = u2->u(x, y);

                x[i] += DX;
                const double fdxdp = u2->u(x, y);

                x[i] = origx - DX;
                vp[j] = origvp;
                u2->setVP(vp);
                const double fmdx = u2->u(x, y);

                vp[j] = origvp + DX;
                u2->setVP(vp);
                const double fmdxdp = u2->u(x,y);

                const double numderiv = (fdxdp - 2.*fdp + fmdxdp - fdx + 2.*f - fmdx)/(DX*DX*DX);

                // cout << "u2->getD2VD1(" << i << ", " << j << ") = " << u2->getD2VD1(i, j) << endl;
                // cout << "numderiv = " << numderiv << endl << endl;
                assert( fabs( (u2->getD2VD1(i, j)-numderiv)/numderiv ) < TINY );

                x[i] = origx;
                vp[j] = origvp;
                u2->setVP(vp);
            }
        }
        for (int i=0; i<NSPACEDIM; ++i){
            for (int j=0; j<poly_u2->getNVP(); ++j){
                const double origx = x[i];
                const double origvp = poly_vp[j];

                x[i] = origx + DX;
                const double fdx = poly_u2->u(x, y);

                x[i] = origx;
                poly_vp[j] += DX;
                poly_u2->setVP(poly_vp);
                const double fdp = poly_u2->u(x, y);

                x[i] += DX;
                const double fdxdp = poly_u2->u(x, y);

                x[i] = origx - DX;
                poly_vp[j] = origvp;
                poly_u2->setVP(poly_vp);
                const double fmdx = poly_u2->u(x, y);

                poly_vp[j] = origvp + DX;
                poly_u2->setVP(poly_vp);
                const double fmdxdp = poly_u2->u(x,y);

                const double numderiv = (fdxdp - 2.*fdp + fmdxdp - fdx + 2.*pf - fmdx)/(DX*DX*DX);

                // cout << "poly_u2->getD2VD1(" << i << ", " << j << ") = " << poly_u2->getD2VD1(i, j) << endl;
                // cout << "numderiv = " << numderiv << endl << endl;
                assert( fabs( (poly_u2->getD2VD1(i, j)-numderiv)/numderiv ) < TINY );

                x[i] = origx;
                poly_vp[j] = origvp;
                poly_u2->setVP(poly_vp);
            }
        }

        // derivative in respect to y
        for (int i=0; i<NSPACEDIM; ++i){
            for (int j=0; j<u2->getNVP(); ++j){
                const double origy = y[i];
                const double origvp = vp[j];

                y[i] = origy + DX;
                const double fdx = u2->u(x, y);

                y[i] = origy;
                vp[j] += DX;
                u2->setVP(vp);
                const double fdp = u2->u(x, y);

                y[i] += DX;
                const double fdxdp = u2->u(x, y);

                y[i] = origy - DX;
                vp[j] = origvp;
                u2->setVP(vp);
                const double fmdx = u2->u(x, y);

                vp[j] = origvp + DX;
                u2->setVP(vp);
                const double fmdxdp = u2->u(x,y);

                const double numderiv = (fdxdp - 2.*fdp + fmdxdp - fdx + 2.*f - fmdx)/(DX*DX*DX);

                // cout << "u2->getD2VD1(" << i+NSPACEDIM << ", " << j << ") = " << u2->getD2VD1(i+NSPACEDIM, j) << endl;
                // cout << "numderiv = " << numderiv << endl << endl;
                assert( fabs( (u2->getD2VD1(i+NSPACEDIM, j)-numderiv)/numderiv ) < TINY );

                y[i] = origy;
                vp[j] = origvp;
                u2->setVP(vp);
            }
        }
        for (int i=0; i<NSPACEDIM; ++i){
            for (int j=0; j<poly_u2->getNVP(); ++j){
                const double origy = y[i];
                const double origvp = poly_vp[j];

                y[i] = origy + DX;
                const double fdx = poly_u2->u(x, y);

                y[i] = origy;
                poly_vp[j] += DX;
                poly_u2->setVP(poly_vp);
                const double fdp = poly_u2->u(x, y);

                y[i] += DX;
                const double fdxdp = poly_u2->u(x, y);

                y[i] = origy - DX;
                poly_vp[j] = origvp;
                poly_u2->setVP(poly_vp);
                const double fmdx = poly_u2->u(x, y);

                poly_vp[j] = origvp + DX;
                poly_u2->setVP(poly_vp);
                const double fmdxdp = poly_u2->u(x,y);

                const double numderiv = (fdxdp - 2.*fdp + fmdxdp - fdx + 2.*pf - fmdx)/(DX*DX*DX);

                // cout << "poly_u2->getD2VD1(" << i+NSPACEDIM << ", " << j << ") = " << poly_u2->getD2VD1(i+NSPACEDIM, j) << endl;
                // cout << "numderiv = " << numderiv << endl << endl;
                assert( fabs( (poly_u2->getD2VD1(i+NSPACEDIM, j)-numderiv)/numderiv ) < TINY );

                y[i] = origy;
                poly_vp[j] = origvp;
                poly_u2->setVP(poly_vp);
            }
        }
    }

    delete poly_u2;
    delete u2;
    delete em;

    return 0;
}
