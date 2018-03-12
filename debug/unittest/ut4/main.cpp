#include "EuclideanMetric.hpp"
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
    const double DX = 0.001;
    const double TINY = 0.01;
    const int NTEST = 10;

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd1, rd2;
    rgen = std::mt19937_64(rdev());
    rgen.seed(18984687);
    rd1 = uniform_real_distribution<double>(-1., -0.1);
    rd2 = uniform_real_distribution<double>(0.1, 1.);

    // distance metric and two body-pseudopotential
    EuclideanMetric * em = new EuclideanMetric(NSPACEDIM);
    He3u2 * u2 = new He3u2(em);

    // position of two particles
    double * x = new double[NSPACEDIM];
    double * y = new double[NSPACEDIM];

    // variational parameters
    double * vp = new double[u2->getNVP()];
    u2->getVP(vp);



    // do NTEST tests with random x and y
    for (int k=0; k<NTEST; ++k){
        // pick random x and y
        for (int i=0; i<NSPACEDIM; ++i){
            x[i] = rd1(rgen);
            y[i] = rd2(rgen);
        }


        // compute value of the u2
        const double f = u2->u(x, y);


        // compute all the analytical derivative
        u2->computeAllDerivatives(x, y);



        // --- check the first derivatives

        // derivative in respect to x
        for (int i=0; i<NSPACEDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            const double fdx = u2->u(x, y);
            const double numderiv = (fdx - f)/DX;

            // cout << "u2->getD1(" << i << ") = " << u2->getD1(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs(u2->getD1(i)-numderiv)/numderiv < TINY );

            x[i] = origx;
        }

        // derivative in respect to y
        for (int i=0; i<NSPACEDIM; ++i){
            const double origy = y[i];
            y[i] += DX;
            const double fdy = u2->u(x, y);
            const double numderiv = (fdy - f)/DX;

            // cout << "u2->getD1(" << i+NSPACEDIM << ") = " << u2->getD2(i+NSPACEDIM) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs(u2->getD1(i+NSPACEDIM)-numderiv)/numderiv < TINY );

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
            assert( abs( (u2->getD2(i)-numderiv)/numderiv ) < TINY );

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
            assert( abs( (u2->getD2(i+NSPACEDIM)-numderiv)/numderiv ) < TINY );

            y[i] = origy;
        }



        // -- check the first variational derivative

        for (int i=0; i<u2->getNVP(); ++i){
            const double origvp = vp[i];
            vp[i] += DX;
            u2->setVP(vp);
            const double fdp = u2->u(x,y);
            const double numderiv = (fdp - f)/DX;

            // cout << "u2->getVD1(" << i << ") = " << u2->getVD1(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( abs(( u2->getVD1(i)-numderiv)/numderiv ) < TINY );

            vp[i] = origvp;
            u2->setVP(vp);
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
                assert( abs( (u2->getD1VD1(i, j)-numderiv)/numderiv ) < TINY );

                x[i] = origx;
                vp[j] = origvp;
                u2->setVP(vp);
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
                assert( abs( (u2->getD1VD1(i+NSPACEDIM, j)-numderiv)/numderiv ) < TINY );

                y[i] = origy;
                vp[j] = origvp;
                u2->setVP(vp);
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
                assert( abs( (u2->getD2VD1(i, j)-numderiv)/numderiv ) < TINY );

                x[i] = origx;
                vp[j] = origvp;
                u2->setVP(vp);
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
                assert( abs( (u2->getD2VD1(i+NSPACEDIM, j)-numderiv)/numderiv ) < TINY );

                y[i] = origy;
                vp[j] = origvp;
                u2->setVP(vp);
            }
        }


    }




    delete[] vp;
    delete[] x;
    delete[] y;
    delete u2;
    delete em;

    return 0;
}
