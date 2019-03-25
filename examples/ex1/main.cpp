#include <cmath>
#include <iostream>

#include "vmc/Hamiltonian.hpp"
#include "vmc/MPIVMC.hpp"
#include "vmc/VMC.hpp"

#include "../common/ExampleFunctions.hpp" // look here for the used WFs / Hamiltonian


int main()
{
    using namespace std;
    using namespace vmc;

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
    vmc->computeEnergy(E_NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    // Case 2
    cout << "-> psi2: " << endl;
    delete vmc;
    vmc = new VMC(psi2, ham2);
    vmc->computeEnergy(E_NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    // Case 3
    cout << "-> psi3: " << endl;
    delete vmc;
    vmc = new VMC(psi3, ham3);
    vmc->computeEnergy(E_NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    // Case 4
    cout << "-> psi4: " << endl;
    delete vmc;
    vmc = new VMC(psi4, ham4);
    vmc->computeEnergy(E_NMC, energy, d_energy);
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
