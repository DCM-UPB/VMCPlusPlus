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

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd;
    rgen = mt19937_64(rdev());
    rgen.seed(18984687);
    rd = uniform_real_distribution<double>(-0.05, 0.05);

    // define one asymmetric and one (naturally) symmetric Wave Function
    const double ai_nosym[NPART] = {0.5, -0.25, 0.0};
    QuadrExponential1DNPOrbital * phi_nosym = new QuadrExponential1DNPOrbital(NPART, ai_nosym, 1.0);

    const double ai_sym[NPART] = {0.25, 0.25, 0.25};
    QuadrExponential1DNPOrbital * psi_sym = new QuadrExponential1DNPOrbital(NPART, ai_sym, 1.0);

    SymmetrizerWaveFunction * phi_sym = new SymmetrizerWaveFunction(phi_nosym, false);

    // particles position
    double x[NPART*NSPACEDIM];
    x[0] = 0.2; x[1] = -0.5; x[2] = 0.7;

    // --- check get/setVP

    // --- check the sampling function

    // --- check the acceptance function
    // slightly change x

    // compute the new protovalues
    double protov;
    //Psi->samplingFunction(x, &protov);

    // compute and compare the acceptance values
    //const double accPsi = Psi->getAcceptance(&protov, &protov);

    // --- check the derivatives

    // pre-compute all the derivatives analytically
    //Psi->computeAllDerivatives(x);

    delete phi_sym;
    delete phi_nosym;
    delete psi_sym;

    return 0;
}
