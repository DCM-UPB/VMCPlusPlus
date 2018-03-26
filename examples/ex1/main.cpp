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
public:
    QuadrExponential1D1POrbital(const double a, const double b):
    WaveFunction(1 /*num space dimensions*/, 1 /*num particles*/, 1 /*num wf components*/, 2 /*num variational parameters*/, true /*VD1*/, false /*D1VD1*/, false /*D2VD1*/) {
            setVP(0, a);
            setVP(1, b);
        }

    void samplingFunction(const double *x, double *out){
        /*
          Compute the sampling function proto value, used in getAcceptance()
        */
        *out = -2.*(getVP(1)*(x[0]-getVP(0))*(x[0]-getVP(0)));
    }

    double getAcceptance(const double * protoold, const double * protonew){
        /*
          Compute the acceptance probability
        */
        return exp(protonew[0]-protoold[0]);
    }

    void computeAllDerivatives(const double *x){
        _setD1DivByWF(0, -2.*getVP(1)*(x[0]-getVP(0)));
        _setD2DivByWF(0, -2.*getVP(1) + (-2.*getVP(1)*(x[0]-getVP(0)))*(-2.*getVP(1)*(x[0]-getVP(0))));
        if (hasVD1()){
            _setVD1DivByWF(0, 2.*getVP(1)*(x[0]-getVP(0)));
            _setVD1DivByWF(1, -(x[0]-getVP(0))*(x[0]-getVP(0)));
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
public:
    Gaussian1D1POrbital(const double b):
    WaveFunction(1, 1, 1, 1, false, false, false){
        setVP(0, b);
    }

    void samplingFunction(const double *in, double *out)
    {
        *out=-2.*getVP(0)*(*in)*(*in);
    }

    double getAcceptance(const double * protoold, const double * protonew)
    {
        return exp(protonew[0]-protoold[0]);
    }

    void computeAllDerivatives(const double *in){
        _setD1DivByWF(0, -2.*getVP(0)*(*in));
        _setD2DivByWF(0, -2.*getVP(0)+4.*getVP(0)*getVP(0)*(*in)*(*in));
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
