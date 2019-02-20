#include "vmc/EuclideanMetric.hpp"
#include "vmc/TwoBodyPseudoPotential.hpp"
#include "vmc/TwoBodyJastrow.hpp"
#include "vmc/MultiComponentWaveFunction.hpp"

#include <assert.h>
#include <math.h>
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
    const double FOO1 = 37.;
    const double FOO2 = 42.;

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd;
    rgen = mt19937_64(rdev());
    rgen.seed(18984687);
    rd = uniform_real_distribution<double>(-0.05, 0.05);

    // Define 4 Jastrow
    EuclideanMetric * em = new EuclideanMetric(NSPACEDIM);
    PolynomialU2 * u2_1 = new PolynomialU2(em, -0.3, -0.1);;
    PolynomialU2 * u2_2 = new PolynomialU2(em, -0.2, -0.15);
    FlatU2 * u2_3 = new FlatU2(em, 3.);
    He3u2 * u2_4 = new He3u2(em);
    TwoBodyJastrow * J_1 = new TwoBodyJastrow(NPART, u2_1);
    TwoBodyJastrow * J_2 = new TwoBodyJastrow(NPART, u2_2);
    TwoBodyJastrow * J_3 = new TwoBodyJastrow(NPART, u2_3);
    TwoBodyJastrow * J_4 = new TwoBodyJastrow(NPART, u2_4);
    TwoBodyJastrow ** J = new TwoBodyJastrow*[4];
    J[0] = J_1; J[1] = J_2; J[2] = J_3; J[3] = J_4;

    // define Multi Component Wave Function
    MultiComponentWaveFunction * Psi = new MultiComponentWaveFunction(NSPACEDIM, NPART);
    Psi->addWaveFunction(J_1);
    Psi->addWaveFunction(J_2);
    Psi->addWaveFunction(J_3);
    Psi->addWaveFunction(J_4);

    // particles position
    double x[NPART*NSPACEDIM];

    // pick x from a grid
    const double K = 0.5;
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
    double vp[Psi->getNVP()];
    Psi->getVP(vp);



    // --- check get/setVP
    // cout << "Psi->getNVP() = " << Psi->getNVP() << endl;
    assert( Psi->getNVP() == J_1->getNVP()+J_2->getNVP()+J_3->getNVP()+J_4->getNVP());
    int contivp = 0;
    for (int iJ=0; iJ<4; ++iJ){
        double Jvp[J[iJ]->getNVP()];
        J[iJ]->getVP(Jvp);
        for (int ivp=0; ivp<J[iJ]->getNVP(); ++ivp){
            // cout << "vp[" << contivp+ivp << "] = " << vp[contivp+ivp] << "    Jvp[" << ivp << "] = " << Jvp[ivp] << endl;
            assert( vp[contivp+ivp] == Jvp[ivp] );
            const double origvp = Jvp[ivp];

            vp[contivp+ivp] = FOO1;
            Psi->setVP(vp);
            Psi->getVP(vp);
            J[iJ]->getVP(Jvp);
            // cout << "vp[" << contivp+ivp << "] = " << vp[contivp+ivp] << "    Jvp[" << ivp << "] = " << Jvp[ivp] << endl;
            assert( vp[contivp+ivp] == Jvp[ivp] );

            Jvp[ivp] = FOO2;
            J[iJ]->setVP(Jvp);
            J[iJ]->getVP(Jvp);
            Psi->getVP(vp);
            // cout << "vp[" << contivp+ivp << "] = " << vp[contivp+ivp] << "    Jvp[" << ivp << "] = " << Jvp[ivp] << endl;
            assert( vp[contivp+ivp] == Jvp[ivp] );

            vp[contivp+ivp] = origvp;
            Jvp[ivp] = origvp;
            Psi->setVP(vp);
        }
        contivp += J[iJ]->getNVP();
    }



    // --- check the sampling function
    double protov[4];
    Psi->samplingFunction(x, protov);
    double * protovJ[4]; for (int i=0; i<4; ++i) protovJ[i] = new double;
    for (int iJ=0; iJ<4; ++iJ){
        J[iJ]->samplingFunction(x, protovJ[iJ]);
        // cout << "Psi protovalue = " << protov[iJ] << "    J_" << iJ+1 << " protovalue = " << protovJ[iJ][0] << endl;
        assert( protov[iJ] == protovJ[iJ][0] );
    }



    // --- check the acceptance function
    // slightly change x
    for (int i=0; i<NPART; ++i){
        for (int j=0; j<NSPACEDIM; ++j){
            x[i*NSPACEDIM+j] += 0.1 * rd(rgen);
        }
    }
    // compute the new protovalues
    double protovnew[4];
    Psi->samplingFunction(x, protovnew);
    double * protovJnew[4]; for (int i=0; i<4; ++i) protovJnew[i] = new double;
    J_1->samplingFunction(x, protovJnew[0]);
    J_2->samplingFunction(x, protovJnew[1]);
    J_3->samplingFunction(x, protovJnew[2]);
    J_4->samplingFunction(x, protovJnew[3]);

    // compute and compare the acceptance values
    const double accPsi = Psi->getAcceptance(protov, protovnew);
    const double accJ1 = J_1->getAcceptance(protovJ[0], protovJnew[0]);
    const double accJ2 = J_2->getAcceptance(protovJ[1], protovJnew[1]);
    const double accJ3 = J_3->getAcceptance(protovJ[2], protovJnew[2]);
    const double accJ4 = J_4->getAcceptance(protovJ[3], protovJnew[3]);
    // cout << "acceptance values:    Psi = " << accPsi << "    J_1 = " << accJ1 << "    J_2 = " << accJ2 << "    J_3 = " << accJ3 << "    J_4 = " << accJ4 << "    J_1*J_2*J_3*J_4 = " << accJ1*accJ2*accJ3*accJ4 << endl;
    assert( accPsi == accJ1*accJ2*accJ3*accJ4 );



    // --- check the derivatives

    // pre-compute all the derivatives analytically
    Psi->computeAllDerivatives(x);


    // initial wave function
    double f, fdx, fmdx, fdvp, fdxdvp, fmdxdvp;
    double samp[4];
    Psi->samplingFunction(x, samp); f = exp(samp[0]+samp[1]+samp[2]+samp[3]);


    // --- check the first derivatives
    for (int i=0; i<NPART*NSPACEDIM; ++i){
        const double origx = x[i];
        x[i] += DX;
        Psi->samplingFunction(x, samp); fdx = exp(samp[0]+samp[1]+samp[2]+samp[3]);
        const double numderiv = (fdx-f)/(DX*f);

        // cout << "getD1DivByWF(" << i <<") = " << Psi->getD1DivByWF(i) << endl;
        // cout << "numderiv = " << numderiv << endl << endl;
        assert( fabs( (Psi->getD1DivByWF(i) - numderiv)/numderiv) < TINY );

        x[i] = origx;
    }


    // --- check the second derivatives
    for (int i=0; i<NPART*NSPACEDIM; ++i){
        const double origx = x[i];
        x[i] = origx + DX;
        Psi->samplingFunction(x, samp); fdx = exp(samp[0]+samp[1]+samp[2]+samp[3]);
        x[i] = origx - DX;
        Psi->samplingFunction(x, samp); fmdx = exp(samp[0]+samp[1]+samp[2]+samp[3]);
        const double numderiv = (fdx - 2.*f + fmdx) / (DX*DX*f);

        // cout << "getD2DivByWF(" << i << ") = " << Psi->getD2DivByWF(i) << endl;
        // cout << "J_1->getD2DivByWF(" << i << ") = " << J_1->getD2DivByWF(i) << endl;
        // cout << "numderiv = " << numderiv << endl << endl;
        assert( fabs( (Psi->getD2DivByWF(i) - numderiv)/numderiv) < TINY );

        x[i] = origx;
    }


    // -- check the first variational derivative
    for (int i=0; i<Psi->getNVP(); ++i){
        const double origvp = vp[i];
        vp[i] += DX;
        Psi->setVP(vp);
        Psi->samplingFunction(x, samp); fdvp = exp(samp[0]+samp[1]+samp[2]+samp[3]);
        const double numderiv = (fdvp - f)/(DX*f);

        // cout << "getVD1DivByWF(" << i << ") = " << Psi->getVD1DivByWF(i) << endl;
        // cout << "numderiv = " << numderiv << endl << endl;
        assert( fabs( (Psi->getVD1DivByWF(i) - numderiv)/numderiv ) < TINY );

        vp[i] = origvp;
        Psi->setVP(vp);
    }


    // --- check the first cross derivative
    for (int i=0; i<NPART*NSPACEDIM; ++i){
        for (int j=0; j<Psi->getNVP(); ++j){
            const double origx = x[i];
            const double origvp = vp[j];

            x[i] += DX;
            Psi->samplingFunction(x, samp); fdx = exp(samp[0]+samp[1]+samp[2]+samp[3]);

            x[i] = origx;
            vp[j] += DX;
            Psi->setVP(vp);
            Psi->samplingFunction(x, samp); fdvp = exp(samp[0]+samp[1]+samp[2]+samp[3]);

            x[i] += DX;
            Psi->samplingFunction(x, samp); fdxdvp = exp(samp[0]+samp[1]+samp[2]+samp[3]);

            const double numderiv = (fdxdvp - fdx - fdvp + f)/(DX*DX*f);

            // cout << "getD1VD1DivByWF(" << i << ", " << j << ") = " << Psi->getD1VD1DivByWF(i, j) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs( (Psi->getD1VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );

            x[i] = origx;
            vp[j] = origvp;
            Psi->setVP(vp);
        }
    }


    // --- check the second cross derivative
    for (int i=0; i<NPART*NSPACEDIM; ++i){
        for (int j=0; j<Psi->getNVP(); ++j){
            const double origx = x[i];
            const double origvp = vp[j];

            x[i] = origx + DX;
            Psi->samplingFunction(x, samp); fdx = exp(samp[0]+samp[1]+samp[2]+samp[3]);

            vp[j] = origvp + DX;
            Psi->setVP(vp);
            Psi->samplingFunction(x, samp); fdxdvp = exp(samp[0]+samp[1]+samp[2]+samp[3]);

            x[i] = origx;
            Psi->samplingFunction(x, samp); fdvp = exp(samp[0]+samp[1]+samp[2]+samp[3]);

            x[i] = origx - DX;
            Psi->samplingFunction(x, samp); fmdxdvp = exp(samp[0]+samp[1]+samp[2]+samp[3]);

            vp[j] = origvp;
            Psi->setVP(vp);
            Psi->samplingFunction(x, samp); fmdx = exp(samp[0]+samp[1]+samp[2]+samp[3]);

            const double numderiv = (fdxdvp - 2.*fdvp + fmdxdvp - fdx + 2.*f - fmdx)/(DX*DX*DX*f);

            // cout << "getD2VD1DivByWF(" << i << ", " << j << ") = " << Psi->getD2VD1DivByWF(i, j) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs( (Psi->getD2VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );

            x[i] = origx;
            vp[j] = origvp;
            Psi->setVP(vp);
        }
    }


    for(int i=0; i<4; ++i){
        delete protovJnew[i];
    }
    for(int i=0; i<4; ++i){
        delete protovJ[i];
    }
    delete Psi;
    delete[] J;
    delete J_4;
    delete J_3;
    delete J_2;
    delete J_1;
    delete u2_4;
    delete u2_3;
    delete u2_2;
    delete u2_1;
    delete em;

    return 0;
}