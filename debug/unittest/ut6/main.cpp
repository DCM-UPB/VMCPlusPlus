#include "EuclideanMetric.hpp"
#include "TwoBodyPseudoPotential.hpp"
#include "TwoBodyJastrow.hpp"
#include "MultiComponentWaveFunction.hpp"

#include <assert.h>
#include <math.h>
#include <iostream>
#include <random>



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
        return _a*r*r + _b*r*r*r;
    }

    double urD1(const double &r){
        return 2.*_a*r + 3.*_b*r*r;
    }

    double urD2(const double &r){
        return 2.*_a + 6.*_b*r;
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
    const int NPART = 3;
    const double DX = 0.001;
    const double TINY = 0.1;
    const double FOO1 = 37.;
    const double FOO2 = 42.;

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd;
    rgen = mt19937_64(rdev());
    rgen.seed(18984687);
    rd = uniform_real_distribution<double>(-0.3, 0.3);

    // Define 4 Jastrow
    EuclideanMetric * em = new EuclideanMetric(NSPACEDIM);
    PolynomialU2 * u2_1 = new PolynomialU2(em, -0.2, -0.2);
    PolynomialU2 * u2_2 = new PolynomialU2(em, -0.1, -0.2);
    PolynomialU2 * u2_3 = new PolynomialU2(em, -0.1, -0.1);
    PolynomialU2 * u2_4 = new PolynomialU2(em, -0.2, -0.1);
    TwoBodyJastrow * J_1 = new TwoBodyJastrow(NPART, u2_1);
    TwoBodyJastrow * J_2 = new TwoBodyJastrow(NPART, u2_2);
    TwoBodyJastrow * J_3 = new TwoBodyJastrow(NPART, u2_3);
    TwoBodyJastrow * J_4 = new TwoBodyJastrow(NPART, u2_4);
    TwoBodyJastrow ** J = new TwoBodyJastrow *[4];
    J[0] = J_1; J[1] = J_2; J[2] = J_3; J[3] = J_4;

    // define Multi Component Wave Function
    MultiComponentWaveFunction * Psi = new MultiComponentWaveFunction(NSPACEDIM, NPART);
    Psi->addWaveFunction(J_1);
    Psi->addWaveFunction(J_2);
    Psi->addWaveFunction(J_3);
    Psi->addWaveFunction(J_4);

    // particles position
    double * x = new double[NPART*NSPACEDIM];

    // pick x from a grid
    const double K = 0.7;
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
    double * vp = new double[Psi->getNVP()];
    Psi->getVP(vp);



    // --- check get/setVP
    cout << "Psi->getNVP() = " << Psi->getNVP() << endl;
    assert( Psi->getNVP() == 8);
    for (int iJ=0; iJ<4; ++iJ){
        double * Jvp = new double[J[iJ]->getNVP()];
        J[iJ]->getVP(Jvp);
        for (int ivp=0; ivp<J[iJ]->getNVP(); ++ivp){
            cout << "vp[" << iJ*2+ivp << "] = " << vp[iJ*2+ivp] << "    Jvp[" << ivp << "] = " << Jvp[ivp] << endl;
            assert( vp[iJ*2+ivp] == Jvp[ivp] );
            const double origvp = Jvp[ivp];

            vp[iJ*2+ivp] = FOO1;
            Psi->setVP(vp);
            Psi->getVP(vp);
            J[iJ]->getVP(Jvp);
            cout << "vp[" << iJ*2+ivp << "] = " << vp[iJ*2+ivp] << "    Jvp[" << ivp << "] = " << Jvp[ivp] << endl;
            assert( vp[iJ*2+ivp] == Jvp[ivp] );

            Jvp[ivp] = FOO2;
            J[iJ]->setVP(Jvp);
            J[iJ]->getVP(Jvp);
            Psi->getVP(vp);
            cout << "vp[" << iJ*2+ivp << "] = " << vp[iJ*2+ivp] << "    Jvp[" << ivp << "] = " << Jvp[ivp] << endl;
            assert( vp[iJ*2+ivp] == Jvp[ivp] );

            vp[iJ*2+ivp] = origvp;
            Jvp[ivp] = origvp;
        }
        delete[] Jvp;
    }



    // --- check the sampling function
    double * protov = new double[4];
    Psi->samplingFunction(x, protov);
    double * protovJ = new double[0];
    for (int iJ=0; iJ<4; ++iJ){
        J[iJ]->samplingFunction(x, protovJ+iJ);
        cout << "Psi protovalue = " << protov[iJ] << "    J_" << iJ+1 << " protovalue = " << protovJ[iJ] << endl;
        assert( protov[iJ] == protovJ[iJ] );
    }



    // --- check the acceptance function
    // slightly change x
    for (int i=0; i<NPART; ++i){
        for (int j=0; j<NSPACEDIM; ++j){
            x[i*NSPACEDIM+j] += 0.1 * rd(rgen);
        }
    }
    // compute the new protovalues
    double * protovnew = new double[4];
    Psi->samplingFunction(x, protovnew);
    double * protovJnew = new double[4];
    J_1->samplingFunction(x, protovJnew);
    J_2->samplingFunction(x, protovJnew+1);
    J_3->samplingFunction(x, protovJnew+2);
    J_4->samplingFunction(x, protovJnew+3);

    // compute and compare the acceptance values
    const double accPsi = Psi->getAcceptance(protov, protovnew);
    const double accJ1 = J_1->getAcceptance(protovJ, protovJnew);
    const double accJ2 = J_2->getAcceptance(protovJ+1, protovJnew+1);
    const double accJ3 = J_3->getAcceptance(protovJ+2, protovJnew+2);
    const double accJ4 = J_4->getAcceptance(protovJ+3, protovJnew+3);
    cout << "acceptance values:    Psi = " << accPsi << "    J_1 = " << accJ1 << "    J_2 = " << accJ2 << "    J_3 = " << accJ3 << "    J_4 = " << accJ4 << "    J_1*J_2*J_3*J_4 = " << accJ1*accJ2*accJ3*accJ4 << endl;
    assert( accPsi == accJ1*accJ2*accJ3*accJ4 );



    // --- check the derivatives

    // pre-compute all the derivatives analytically
    Psi->computeAllDerivatives(x);


    // // initial wave function
    // double f, fdx, fmdx, fdvp, fdxdvp, fmdxdvp;
    // J->samplingFunction(x, &f); f = exp(f);


    // // --- check the first derivatives
    // for (int i=0; i<NPART*NSPACEDIM; ++i){
    //     const double origx = x[i];
    //     x[i] += DX;
    //     J->samplingFunction(x, &fdx); fdx = exp(fdx);
    //     const double numderiv = (fdx-f)/(DX*f);
    //
    //     // cout << "getD1DivByWF(" << i <<") = " << J->getD1DivByWF(i) << endl;
    //     // cout << "numderiv = " << numderiv << endl << endl;
    //     assert( abs( (J->getD1DivByWF(i) - numderiv)/numderiv) < TINY );
    //
    //     x[i] = origx;
    // }
    //
    //
    // // --- check the second derivatives
    // for (int i=0; i<NPART*NSPACEDIM; ++i){
    //     const double origx = x[i];
    //     x[i] += DX;
    //     J->samplingFunction(x, &fdx); fdx = exp(fdx);
    //     x[i] -= 2.*DX;
    //     J->samplingFunction(x, &fmdx); fmdx = exp(fmdx);
    //     const double numderiv = (fdx - 2.*f + fmdx) / (DX*DX*f);
    //
    //     // cout << "getD2DivByWF(" << i << ") = " << J->getD2DivByWF(i) << endl;
    //     // cout << "numderiv = " << numderiv << endl << endl;
    //     assert( abs( (J->getD2DivByWF(i) - numderiv)/numderiv) < TINY );
    //
    //     x[i] = origx;
    // }
    //
    //
    // // -- check the first variational derivative
    // for (int i=0; i<J->getNVP(); ++i){
    //     const double origvp = vp[i];
    //     vp[i] += DX;
    //     J->setVP(vp);
    //     J->samplingFunction(x, &fdvp); fdvp = exp(fdvp);
    //     const double numderiv = (fdvp - f)/(DX*f);
    //
    //     // cout << "getVD1DivByWF(" << i << ") = " << J->getVD1DivByWF(i) << endl;
    //     // cout << "numderiv = " << numderiv << endl << endl;
    //     assert( abs( (J->getVD1DivByWF(i) - numderiv)/numderiv ) < TINY );
    //
    //     vp[i] = origvp;
    //     J->setVP(vp);
    // }
    //
    //
    // // --- check the first cross derivative
    // for (int i=0; i<NPART*NSPACEDIM; ++i){
    //     for (int j=0; j<J->getNVP(); ++j){
    //         const double origx = x[i];
    //         const double origvp = vp[j];
    //
    //         x[i] += DX;
    //         J->samplingFunction(x, &fdx); fdx = exp(fdx);
    //
    //         x[i] = origx;
    //         vp[j] += DX;
    //         J->setVP(vp);
    //         J->samplingFunction(x, &fdvp); fdvp = exp(fdvp);
    //
    //         x[i] += DX;
    //         J->samplingFunction(x, &fdxdvp); fdxdvp = exp(fdxdvp);
    //
    //         const double numderiv = (fdxdvp - fdx - fdvp + f)/(DX*DX*f);
    //
    //         // cout << "getD1VD1DivByWF(" << i << ", " << j << ") = " << J->getD1VD1DivByWF(i, j) << endl;
    //         // cout << "numderiv = " << numderiv << endl << endl;
    //         assert( abs( (J->getD1VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );
    //
    //         x[i] = origx;
    //         vp[j] = origvp;
    //         J->setVP(vp);
    //     }
    // }
    //
    //
    // // --- check the second cross derivative
    // for (int i=0; i<NPART*NSPACEDIM; ++i){
    //     for (int j=0; j<J->getNVP(); ++j){
    //         const double origx = x[i];
    //         const double origvp = vp[j];
    //
    //         x[i] = origx + DX;
    //         J->samplingFunction(x, &fdx); fdx = exp(fdx);
    //
    //         vp[j] = origvp + DX;
    //         J->setVP(vp);
    //         J->samplingFunction(x, &fdxdvp); fdxdvp = exp(fdxdvp);
    //
    //         x[i] = origx;
    //         J->samplingFunction(x, &fdvp); fdvp = exp(fdvp);
    //
    //         x[i] = origx - DX;
    //         J->samplingFunction(x, &fmdxdvp); fmdxdvp = exp(fmdxdvp);
    //
    //         vp[j] = origvp;
    //         J->setVP(vp);
    //         J->samplingFunction(x, &fmdx); fmdx = exp(fmdx);
    //
    //         const double numderiv = (fdxdvp - 2.*fdvp + fmdxdvp - fdx + 2.*f - fmdx)/(DX*DX*DX*f);
    //
    //         // cout << "getD2VD1DivByWF(" << i << ", " << j << ") = " << J->getD2VD1DivByWF(i, j) << endl;
    //         // cout << "numderiv = " << numderiv << endl << endl;
    //         assert( abs( (J->getD2VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );
    //
    //         x[i] = origx;
    //         vp[j] = origvp;
    //         J->setVP(vp);
    //     }
    // }




    delete[] protovJnew;
    delete[] protovnew;
    delete[] protovJ;
    delete[] protov;
    delete[] vp;
    delete[] x;
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
