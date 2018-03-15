#include <iostream>
#include <cmath>
#include <stdexcept>

#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "VMC.hpp"
#include "ConjGrad.hpp"



/*
  Hamiltonian describing a 1-particle harmonic oscillator:
  H  =  p^2 / 2m  +  1/2 * w^2 * x^2
*/
class HarmonicOscillator1D1P: public Hamiltonian{

protected:
    double _w;

public:
    HarmonicOscillator1D1P(const double w, WaveFunction * wf):
        Hamiltonian(1 /*num space dimensions*/, 1 /*num particles*/, wf) {_w=w;}
    virtual ~HarmonicOscillator1D1P(){}

    // potential energy
    double localPotentialEnergy(const double *r)
    {
        return (0.5*_w*_w*(*r)*(*r));
    }
};



/*
  Trial Wave Function for 1 particle in 1 dimension, that uses two variational parameters: a and b.
  Psi  =  exp( -b * (x-a)^2 )
  Notice that the corresponding probability density (sampling function) is Psi^2.
*/
class QuadrExponential1D1POrbital: public WaveFunction{
protected:
    double _a, _b;

public:
    QuadrExponential1D1POrbital(const double a, const double b):
    WaveFunction(1 /*num space dimensions*/, 1 /*num particles*/, 1 /*num wf components*/, 2 /*num variational parameters*/, false /*VD1*/, false /*D1VD1*/, false /*D2VD1*/) {
            _a=a; _b=b;
        }

    void setVP(const double *in){
        _a=in[0];
        _b=in[1];
    }

    void getVP(double *out){
        out[0]=_a;
        out[1]=_b;
    }

    void samplingFunction(const double *x, double *out){
        /*
          Compute the sampling function proto value, used in getAcceptance()
        */
        *out = -2.*(_b*(x[0]-_a)*(x[0]-_a));
    }

    double getAcceptance(const double * protoold, const double * protonew){
        /*
          Compute the acceptance probability
        */
        return exp(protonew[0]-protoold[0]);
    }

    void computeAllDerivatives(const double *x){
        _setD1DivByWF(0, -2.*_b*(x[0]-_a));
        _setD2DivByWF(0, -2.*_b + (-2.*_b*(x[0]-_a))*(-2.*_b*(x[0]-_a)));
        if (hasVD1()){
            _setVD1DivByWF(0, 2.*_b*(x[0]-_a));
            _setVD1DivByWF(1, -(x[0]-_a)*(x[0]-_a));
        }
    }

};



/*
  Trial Wave Function for 1 particle in 1 dimension, that uses one variational parameters: b.
  Psi  =  exp( -b * x^2 )
  Notice that the corresponding probability density (sampling function) is Psi^2.
*/
class Gaussian1D1POrbital: public WaveFunction
{
protected:
    double _b;

public:
    Gaussian1D1POrbital(const double b):
    WaveFunction(1, 1, 1, 1, false, false, false){
        _b=b;
    }

    void setVP(const double *in)
    {
        _b=*in;
        //if (_b<0.01) _b=0.01;
        using namespace std;
        //cout << "change b! " << _b << endl;
    }
    void getVP(double *out)
    {
        *out=_b;
    }

    void samplingFunction(const double *in, double *out)
    {
        *out=-2.*_b*(*in)*(*in);
    }

    double getAcceptance(const double * protoold, const double * protonew)
    {
        return exp(protonew[0]-protoold[0]);
    }

    void computeAllDerivatives(const double *in){
        _setD1DivByWF(0, -2.*_b*(*in));
        _setD2DivByWF(0, -2.*_b+4.*_b*_b*(*in)*(*in));
        if (hasVD1()){
            _setVD1DivByWF(0, (-(*in)*(*in)));
        }
    }
};



int main(){
    using namespace std;

    // Declare some trial wave functions
    Gaussian1D1POrbital * psi1 = new Gaussian1D1POrbital(1.2);
    Gaussian1D1POrbital * psi2 = new Gaussian1D1POrbital(0.5);
    Gaussian1D1POrbital * psi3 = new Gaussian1D1POrbital(1.0);
    QuadrExponential1D1POrbital * psi4 = new QuadrExponential1D1POrbital(-0.5, 1.0);

    // Declare an Hamiltonian for each wave function (keep in mind that the kinetic energy is strictly bound to it)
    // We use the harmonic oscillator with w=1
    HarmonicOscillator1D1P * ham1 = new HarmonicOscillator1D1P(1., psi1);
    HarmonicOscillator1D1P * ham2 = new HarmonicOscillator1D1P(1., psi2);
    HarmonicOscillator1D1P * ham3 = new HarmonicOscillator1D1P(1., psi3);
    HarmonicOscillator1D1P * ham4 = new HarmonicOscillator1D1P(1., psi4);



    cout << endl << " - - - EVALUATION OF ENERGY - - - " << endl << endl;

    VMC * vmc; // VMC object we will resuse
    const long E_NMC = 100000l; // MC samplings to use for computing the energy
    double * energy = new double[4]; // energy
    double * d_energy = new double[4]; // energy error bar

    // Case 1
    cout << "-> psi1: " << endl;
    vmc = new VMC(psi1, ham1);
    vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    // Case 2
    cout << "-> psi2: " << endl;
    delete vmc;
    vmc = new VMC(psi2, ham2);
    vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    // Case 3
    cout << "-> psi3: " << endl;
    delete vmc;
    vmc = new VMC(psi3, ham3);
    vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    // Case 4
    cout << "-> psi4: " << endl;
    delete vmc;
    vmc = new VMC(psi4, ham4);
    vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;


    delete[] d_energy;
    delete[] energy;
    delete vmc;


    delete ham4;
    delete ham3;
    delete ham2;
    delete ham1;

    delete psi4;
    delete psi3;
    delete psi2;
    delete psi1;



    return 0;
}
