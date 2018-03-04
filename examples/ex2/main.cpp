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
    WaveFunction(1 /*num space dimensions*/, 1 /*num particles*/, 1 /*num wf components*/, 2 /*num variational parameters*/, true /*VD1*/, false /*D1VD1*/, false /*D2VD1*/) {
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

    double getAcceptance(){
        /*
          Compute the acceptance probability
        */
        return exp(getProtoNew(0)-getProtoOld(0));
    }

    void computeAllDerivatives(const double *x){
        setD1DivByWF(0, -2.*_b*(x[0]-_a));
        setD2DivByWF(0, -2.*_b + (-2.*_b*(x[0]-_a))*(-2.*_b*(x[0]-_a)));
        if (hasVD1()){
            setVD1DivByWF(0, 2.*_b*(x[0]-_a));
            setVD1DivByWF(1, -(x[0]-_a)*(x[0]-_a));
        }
    }

};




int main(){
    using namespace std;

    // Declare some trial wave functions
    QuadrExponential1D1POrbital * psi = new QuadrExponential1D1POrbital(-0.5, 1.0);

    // Declare an Hamiltonian
    // We use the harmonic oscillator with w=1 and w=2
    const double w1 = 1.;
    HarmonicOscillator1D1P * ham1 = new HarmonicOscillator1D1P(w1, psi);
    const double w2 = 2.;
    HarmonicOscillator1D1P * ham2 = new HarmonicOscillator1D1P(w2, psi);



    cout << endl << " - - - WAVE FUNCTION OPTIMIZATION - - - " << endl << endl;

    VMC * vmc; // VMC object we will resuse
    const long E_NMC = 4000l; // MC samplings to use for computing the energy
    const long G_NMC = 10000l; // MC samplings to use for computing the energy gradient
    double * energy = new double[4]; // energy
    double * d_energy = new double[4]; // energy error bar
    double * vp = new double[psi->getNVP()];


    // Case 1
    cout << "-> ham1:    w = " << w1 << endl << endl;
    vmc = new VMC(psi, ham1);

    cout << "   Initial Wave Function parameters:" << endl;
    psi->getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Starting energy:" << endl;
    vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    cout << "   Optimization . . ." << endl;
    vmc->conjugateGradientOptimization(E_NMC, G_NMC);
    cout << "   . . . Done!" << endl << endl;

    cout << "   Optimized Wave Function parameters:" << endl;
    psi->getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Optimized energy:" << endl;
    vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl << endl;



    // Case 2
    cout << "-> ham2:    w = " << w2 << endl << endl;
    delete vmc;
    vmc = new VMC(psi, ham2);

    cout << "   Initial Wave Function parameters:" << endl;
    psi->getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Starting energy:" << endl;
    vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    cout << "   Optimization . . ." << endl;
    vmc->conjugateGradientOptimization(E_NMC, G_NMC);
    cout << "   . . . Done!" << endl << endl;

    cout << "   Optimized Wave Function parameters:" << endl;
    psi->getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Optimized energy:" << endl;
    vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;




    delete[] vp;
    delete[] d_energy;
    delete[] energy;
    delete vmc;

    delete ham2;
    delete ham1;
    delete psi;



    return 0;
}
