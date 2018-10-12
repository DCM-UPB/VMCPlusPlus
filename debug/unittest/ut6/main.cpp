#include "SymmetrizerWaveFunction.hpp"

#include <assert.h>
#include <math.h>
#include <iostream>
#include <random>


/*
  Trial Wave Function for N particles in 1 dimension, that uses parameters ai and b, with fixed ai and variational b.
  Psi  =  exp( -b * sum((xi-ai)^2) )
*/
class QuadrExponential1DNPOrbital: public WaveFunction{
protected:
    double * _a;
    double _b;

public:
    QuadrExponential1DNPOrbital(const int &npart, const double * a, const double &b):
    WaveFunction(1, npart, 1, 1, true, true, true), _b(b)
    {
        _a=new double[npart];
        for (int i=0; i<npart; ++i) _a[i] = a[i];
    }

    ~QuadrExponential1DNPOrbital()
    {
        delete [] _a;
    }

    void setVP(const double *in){
        _b=in[0];
    }

    void getVP(double *out){
        out[0]=_b;
    }

    void samplingFunction(const double *x, double *out)
    {
        *out = 0.;
        for (int i=0; i<_npart; ++i) {
            *out += (x[i]-_a[i])*(x[i]-_a[i]);
        }
        *out *= -2.*_b;
    }

    double getAcceptance(const double * protoold, const double * protonew)
    {
        return exp(protonew[0]-protoold[0]);
    }

    void computeAllDerivatives(const double *x)
    {
        if (hasVD1()){
            _setVD1DivByWF(0, 0);
        }
        for (int i=0; i<_npart; ++i) {
            _setD1DivByWF(i, -2.*_b*(x[i]-_a[i]));
            _setD2DivByWF(i, -2.*_b + (-2.*_b*(x[i]-_a[i]))*(-2.*_b*(x[i]-_a[i])));
            if (hasVD1()){
                _setVD1DivByWF(0, getVD1DivByWF(0) - (x[i]-_a[i])*(x[i]-_a[i]));
            }
            if (hasD1VD1()){
                _setD1VD1DivByWF(i, 0, -2.*(x[i]-_a[i]) + 2.*_b*(x[i]-_a[i])*(x[i]-_a[i])*(x[i]-_a[i]));
            }
            if (hasD2VD1()){
                _setD2VD1DivByWF(i, 0, -2. + 6.*_b*(x[i]-_a[i])*(x[i]-_a[i]) - 4.*_b*_b*(x[i]-_a[i])*(x[i]-_a[i])*(x[i]-_a[i])*(x[i]-_a[i]));
            }
        }
    }

    double computeWFValue(const double * protovalues)
    {
        return exp(0.5*protovalues[0]);
    }
};


int main(){
    using namespace std;

    // constants
    const int NSPACEDIM = 1;
    const int NPART = 3;
    const double DX = 0.0001;
    const double TINY = 0.01;
    const double SUPERTINY = 0.0000001;

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd;
    rgen = mt19937_64(rdev());
    rgen.seed(18984687);
    rd = uniform_real_distribution<double>(-0.05, 0.05);

    // define one not symmetric and one (naturally) symmetric Wave Function
    const double ai_nosym[NPART] = {0.5, -0.25, 0.0};
    QuadrExponential1DNPOrbital * phi_nosym = new QuadrExponential1DNPOrbital(NPART, ai_nosym, 1.0);

    const double ai_sym[NPART] = {0.25, 0.25, 0.25};
    QuadrExponential1DNPOrbital * psi_sym = new QuadrExponential1DNPOrbital(NPART, ai_sym, 1.0);

    SymmetrizerWaveFunction * phi_sym = new SymmetrizerWaveFunction(phi_nosym, false);
    SymmetrizerWaveFunction * phi_asym = new SymmetrizerWaveFunction(phi_nosym, true); // antisymmetric version

    // particles position and all permutations
    double ** xp = new double * [6];
    for (int i=0; i<6; ++i) xp[i] = new double[NPART*NSPACEDIM];
    xp[0][0] = 0.2; xp[0][1] = -0.5; xp[0][2] = 0.7;
    xp[1][0] = xp[0][0]; xp[1][1] = xp[0][2]; xp[1][2] = xp[0][1];
    xp[2][0] = xp[0][1]; xp[2][1] = xp[0][0]; xp[2][2] = xp[0][2];
    xp[3][0] = xp[0][1]; xp[3][1] = xp[0][2]; xp[3][2] = xp[0][0];
    xp[4][0] = xp[0][2]; xp[4][1] = xp[0][0]; xp[4][2] = xp[0][1];
    xp[5][0] = xp[0][2]; xp[5][1] = xp[0][1]; xp[5][2] = xp[0][0];

    bool isOdd[6] = {false, true, true, false, false, true};

    // compute the new protovalues
    double protov_nosym[2], protov_sym[2], protov_asym[2];
    
    phi_nosym->samplingFunction(xp[0], protov_nosym);
    phi_sym->samplingFunction(xp[0], protov_sym);
    phi_asym->samplingFunction(xp[0], protov_asym);
    for (int i=0; i<6; ++i) {
        size_t offset = (i==0 ? 0 : 1);

        phi_nosym->samplingFunction(xp[i], protov_nosym + offset);
        phi_sym->samplingFunction(xp[i], protov_sym + offset);
        phi_asym->samplingFunction(xp[i], protov_asym + offset);

        /*
        cout << "Perm " << i << ": ";
        cout << " phi_nosym " << phi_nosym->computeWFValue(protov_nosym + offset);
        cout << " phi_sym " << phi_sym->computeWFValue(protov_sym + offset);
        cout << " phi_asym " << phi_asym->computeWFValue(protov_asym + offset) << endl;
        */

        if (i>0) {
            assert(protov_nosym[0] != protov_nosym[1]);
            if (isOdd[i]) assert(abs(protov_asym[0] + protov_asym[1]) < SUPERTINY);
            else assert(abs(protov_asym[0] - protov_asym[1]) < SUPERTINY);
            assert(abs(protov_sym[0] - protov_sym[1]) < SUPERTINY);
        }
    }


    // pre-compute all the derivatives analytically
    phi_sym->computeAllDerivatives(xp[0]);
    phi_asym->computeAllDerivatives(xp[0]);

/*
    // initial wave function
    double f_sym, fdx_sym, fmdx_sym, fdvp_sym, fdxdvp_sym, fmdxdvp_sym;
    double * samp = new double[4];
    Psi->samplingFunction(x, samp); f = exp(samp[0]+samp[1]+samp[2]+samp[3]);


    // --- check the first derivatives
    for (int i=0; i<NPART*NSPACEDIM; ++i){
        const double origx = x[i];
        x[i] += DX;
        Psi->samplingFunction(x, samp); fdx = exp(samp[0]+samp[1]+samp[2]+samp[3]);
        const double numderiv = (fdx-f)/(DX*f);

        // cout << "getD1DivByWF(" << i <<") = " << Psi->getD1DivByWF(i) << endl;
        // cout << "numderiv = " << numderiv << endl << endl;
        assert( abs( (Psi->getD1DivByWF(i) - numderiv)/numderiv) < TINY );

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
        assert( abs( (Psi->getD2DivByWF(i) - numderiv)/numderiv) < TINY );

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
        assert( abs( (Psi->getVD1DivByWF(i) - numderiv)/numderiv ) < TINY );

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
            assert( abs( (Psi->getD1VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );

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
            assert( abs( (Psi->getD2VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );

            x[i] = origx;
            vp[j] = origvp;
            Psi->setVP(vp);
        }
    }

*/
    for (int i=0; i<6; ++i) delete [] xp[i];
    delete [] xp;

    delete phi_sym;
    delete phi_asym;

    delete phi_nosym;
    delete psi_sym;

    return 0;
}
