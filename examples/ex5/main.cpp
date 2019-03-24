#include <cmath>
#include <iostream>
#include <stdexcept>

#include "vmc/Hamiltonian.hpp"
#include "vmc/MPIVMC.hpp"
#include "vmc/VMC.hpp"
#include "vmc/WaveFunction.hpp"

#include "../common/ExampleFunctions.hpp"

int main(){
    using namespace std;

    MPIVMC::Init(); // make this usable with a MPI-compiled library

    // Declare some trial wave functions
    auto * psi = new QuadrExponential1D1POrbital(-0.5, 1.0);

    // Declare an Hamiltonian
    // We use the harmonic oscillator with w=1 and w=2
    const double w1 = 1.;
    auto * ham1 = new HarmonicOscillator1D1P(w1, psi);
    const double w2 = 2.;
    auto * ham2 = new HarmonicOscillator1D1P(w2, psi);



    cout << endl << " - - - WAVE FUNCTION OPTIMIZATION - - - " << endl << endl;

    VMC * vmc; // VMC object we will resuse
    const int E_NMC = 10000l; // MC samplings to use for computing the energy
    double energy[4]; // energy
    double d_energy[4]; // energy error bar
    double vp[psi->getNVP()];


    // Case 1
    cout << "-> ham1:    w = " << w1 << endl << endl;
    vmc = new VMC(psi, ham1);

    cout << "   Initial Wave Function parameters:" << endl;
    psi->getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Starting energy:" << endl;
    vmc->computeEnergy(E_NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    cout << "   Optimization . . ." << endl;

    // settings for better performance
    vmc->getMCI()->setNfindMRT2Iterations(10);
    vmc->getMCI()->setNdecorrelationSteps(1000);

    vmc->nmsimplexOptimization(E_NMC, 1.0 /* weight for energy cost*/, 0.1 /* weight for energy error cost*/, 0.005 /* weight for regularization cost */, 0.1 /* initial simplex size */);
    cout << "   . . . Done!" << endl << endl;

    cout << "   Optimized Wave Function parameters:" << endl;
    psi->getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Optimized energy:" << endl;
    vmc->computeEnergy(E_NMC, energy, d_energy);
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
    vmc->computeEnergy(E_NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    cout << "   Optimization . . ." << endl;

    // settings for better performance
    vmc->getMCI()->setNfindMRT2Iterations(10);
    vmc->getMCI()->setNdecorrelationSteps(1000);

    vmc->nmsimplexOptimization(E_NMC, 1.0 /* weight for energy cost*/, 0.1 /* weight for energy error cost*/, 0.005 /* weight for regularization cost */, 0.1 /* initial simplex size */);
    cout << "   . . . Done!" << endl << endl;

    cout << "   Optimized Wave Function parameters:" << endl;
    psi->getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Optimized energy:" << endl;
    vmc->computeEnergy(E_NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;



    delete vmc;

    delete ham2;
    delete ham1;
    delete psi;


    MPIVMC::Finalize();

    return 0;
}
