#include <iostream>
#include <cmath>
#include <stdexcept>

#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "VMC.hpp"
#include "ConjGrad.hpp"

#include "MPIVMC.hpp"

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

    int myrank = MPIVMC::Init();

    // Declare some trial wave functions
    Gaussian1D1POrbital * psi2 = new Gaussian1D1POrbital(0.5);

    // Declare an Hamiltonian for each wave function (keep in mind that the kinetic energy is strictly bound to it)
    // We use the harmonic oscillator with w=1
    HarmonicOscillator1D1P * ham2 = new HarmonicOscillator1D1P(1., psi2);


    cout << endl << " - - - EVALUATION OF ENERGY - - - " << endl << endl;

    VMC * vmc; // VMC object we will resuse
    const long E_NMC = 10000l; // MC samplings to use for computing the energy
    double * energy = new double[4]; // energy
    double * d_energy = new double[4]; // energy error bar

    // Case 2
    cout << "-> psi2: " << endl;
    delete vmc;
    vmc = new VMC(psi2, ham2);

    auto obsfile = "obsfile";
    auto wlkfile = "wlkfile";
    vmc->getMCI()->storeObservablesOnFile(obsfile, 1);
    vmc->getMCI()->storeWalkerPositionsOnFile(wlkfile, 1);

    vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    delete[] d_energy;
    delete[] energy;
    delete vmc;


    delete ham2;

    delete psi2;

    MPIVMC::Finalize();

    return 0;
}
