#include <cmath>
#include <iostream>
#include <stdexcept>

#include "nfm/ConjGrad.hpp"
#include "vmc/Hamiltonian.hpp"
#include "vmc/MPIVMC.hpp"
#include "vmc/VMC.hpp"
#include "vmc/WaveFunction.hpp"


/*
  Hamiltonian describing a 1-particle harmonic oscillator:
  H  =  p^2 / 2m  +  1/2 * w^2 * x^2
*/
class HarmonicOscillator1D1P: public Hamiltonian
{
protected:
    double _w;

public:
    HarmonicOscillator1D1P(const double w, WaveFunction * wf):
        Hamiltonian(1 /*num space dimensions*/, 1 /*num particles*/, wf) {_w=w;}
    ~HarmonicOscillator1D1P() override= default;

    // potential energy
    double localPotentialEnergy(const double *r) override
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

    void setVP(const double *in) override{
        _a=in[0];
        _b=in[1];
    }

    void getVP(double *out) override{
        out[0]=_a;
        out[1]=_b;
    }

    void protoFunction(const double *x, double *out) override{
        /*
          Compute the sampling function proto value, used in acceptanceFunction()
        */
        *out = -2.*(_b*(x[0]-_a)*(x[0]-_a));
    }

    double acceptanceFunction(const double * protoold, const double * protonew) override{
        /*
          Compute the acceptance probability
        */
        return exp(protonew[0]-protoold[0]);
    }

    void computeAllDerivatives(const double *x) override{
        _setD1DivByWF(0, -2.*_b*(x[0]-_a));
        _setD2DivByWF(0, -2.*_b + (-2.*_b*(x[0]-_a))*(-2.*_b*(x[0]-_a)));
        if (hasVD1()){
            _setVD1DivByWF(0, 2.*_b*(x[0]-_a));
            _setVD1DivByWF(1, -(x[0]-_a)*(x[0]-_a));
        }
    }

    double computeWFValue(const double * protovalues) override
    {
        return exp(0.5*protovalues[0]);
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
    explicit Gaussian1D1POrbital(const double b):
    WaveFunction(1, 1, 1, 1, false, false, false){
        _b=b;
    }

    void setVP(const double *in) override
    {
        _b=*in;
        //if (_b<0.01) _b=0.01;
        using namespace std;
        //cout << "change b! " << _b << endl;
    }
    void getVP(double *out) override
    {
        *out=_b;
    }

    void protoFunction(const double *in, double *out) override
    {
        *out=-2.*_b*(*in)*(*in);
    }

    double acceptanceFunction(const double * protoold, const double * protonew) override
    {
        return exp(protonew[0]-protoold[0]);
    }

    void computeAllDerivatives(const double *in) override{
        _setD1DivByWF(0, -2.*_b*(*in));
        _setD2DivByWF(0, -2.*_b+4.*_b*_b*(*in)*(*in));
        if (hasVD1()){
            _setVD1DivByWF(0, (-(*in)*(*in)));
        }
    }

    double computeWFValue(const double * protovalues) override
    {
        return exp(0.5*protovalues[0]);
    }
};



int main(){
    using namespace std;

    MPIVMC::Init(); // make this usable with a MPI-compiled library

    // Declare some trial wave functions
    auto * psi1 = new Gaussian1D1POrbital(1.2);
    auto * psi2 = new Gaussian1D1POrbital(0.5);
    auto * psi3 = new Gaussian1D1POrbital(1.0);
    auto * psi4 = new QuadrExponential1D1POrbital(-0.5, 1.0);

    // Declare an Hamiltonian for each wave function (keep in mind that the kinetic energy is strictly bound to it)
    // We use the harmonic oscillator with w=1
    auto * ham1 = new HarmonicOscillator1D1P(1., psi1);
    auto * ham2 = new HarmonicOscillator1D1P(1., psi2);
    auto * ham3 = new HarmonicOscillator1D1P(1., psi3);
    auto * ham4 = new HarmonicOscillator1D1P(1., psi4);



    cout << endl << " - - - EVALUATION OF ENERGY - - - " << endl << endl;

    VMC * vmc; // VMC object we will resuse
    const int E_NMC = 100000l; // MC samplings to use for computing the energy
    double energy[4]; // energy
    double d_energy[4]; // energy error bar

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


    delete vmc;


    delete ham4;
    delete ham3;
    delete ham2;
    delete ham1;

    delete psi4;
    delete psi3;
    delete psi2;
    delete psi1;

    MPIVMC::Finalize();

    return 0;
}
