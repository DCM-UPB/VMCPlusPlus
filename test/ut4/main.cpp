#include "vmc/EuclideanMetric.hpp"
#include "vmc/TwoBodyJastrow.hpp"
#include "vmc/TwoBodyPseudoPotential.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>

#include "TestVMCFunctions.hpp"


int main(){
    using namespace std;

    // constants
    const int NSPACEDIM = 3;
    const int NPART = 3;
    const double DX = 0.0001;
    const double TINY = 0.01;

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd;
    rgen = mt19937_64(rdev());
    rgen.seed(18984687);
    rd = uniform_real_distribution<double>(-0.1, 0.1);

    // distance metric
    auto * em = new EuclideanMetric(NSPACEDIM);

    // make the tests using two different pseudo-potentials
    for (int iu2=0; iu2<2; ++iu2) {
        // declare the pseudo-potential and the Jastrow
        TwoBodyPseudoPotential * u2;
        TwoBodyJastrow * J;
        if (iu2 == 0){
            u2 = new PolynomialU2(em, -0.3, -0.1);
        } else if (iu2==1){
            u2 = new He3u2(em);
        }
        J = new TwoBodyJastrow(NPART, u2);


        // particles position
        double x[NPART*NSPACEDIM];

        // pick x from a grid
        const double K = 0.75;
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
        double vp[J->getNVP()];
        J->getVP(vp);


        // compute all the derivatives analytically
        J->computeAllDerivatives(x);


        // initial wave function
        double f, fdx, fmdx, fdp, fdxdp, fmdxdp;
        J->samplingFunction(x, &f); f = exp(f);
        // cout << "f = " << f << endl;


        // --- check the first derivatives
        for (int i=0; i<NPART*NSPACEDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            J->samplingFunction(x, &fdx); fdx = exp(fdx);
            const double numderiv = (fdx-f)/(DX*f);

            // cout << "getD1DivByWF(" << i <<") = " << J->getD1DivByWF(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs( (J->getD1DivByWF(i) - numderiv)/numderiv) < TINY );

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

            // cout << "getD2DivByWF(" << i << ") = " << J->getD2DivByWF(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs( (J->getD2DivByWF(i) - numderiv)/numderiv) < TINY );

            x[i] = origx;
        }


        // -- check the first variational derivative
        for (int i=0; i<J->getNVP(); ++i){
            const double origvp = vp[i];
            vp[i] += DX;
            J->setVP(vp);
            J->samplingFunction(x, &fdp); fdp = exp(fdp);
            const double numderiv = (fdp - f)/(DX*f);

            // cout << "getVD1DivByWF(" << i << ") = " << J->getVD1DivByWF(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs( (J->getVD1DivByWF(i) - numderiv)/numderiv ) < TINY );

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
                J->samplingFunction(x, &fdp); fdp = exp(fdp);

                x[i] += DX;
                J->samplingFunction(x, &fdxdp); fdxdp = exp(fdxdp);

                const double numderiv = (fdxdp - fdx - fdp + f)/(DX*DX*f);

                // cout << "getD1VD1DivByWF(" << i << ", " << j << ") = " << J->getD1VD1DivByWF(i, j) << endl;
                // cout << "numderiv = " << numderiv << endl << endl;
                assert( fabs( (J->getD1VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );

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
                J->samplingFunction(x, &fdxdp); fdxdp = exp(fdxdp);

                x[i] = origx;
                J->samplingFunction(x, &fdp); fdp = exp(fdp);

                x[i] = origx - DX;
                J->samplingFunction(x, &fmdxdp); fmdxdp = exp(fmdxdp);

                vp[j] = origvp;
                J->setVP(vp);
                J->samplingFunction(x, &fmdx); fmdx = exp(fmdx);

                const double numderiv = (fdxdp - 2.*fdp + fmdxdp - fdx + 2.*f - fmdx)/(DX*DX*DX*f);

                // cout << "getD2VD1DivByWF(" << i << ", " << j << ") = " << J->getD2VD1DivByWF(i, j) << endl;
                // cout << "numderiv = " << numderiv << endl << endl;
                assert( fabs( (J->getD2VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );

                x[i] = origx;
                vp[j] = origvp;
                J->setVP(vp);
            }
        }

        delete J;
        delete u2;
    }

    delete em;

    return 0;
}
