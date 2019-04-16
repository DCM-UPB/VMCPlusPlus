#include <cmath>
#include <iostream>
#include <memory>

#include "vmc/MPIVMC.hpp"
#include "vmc/VMC.hpp"

#include "../common/ExampleFunctions.hpp" // look here for the used WFs / Hamiltonian


int main()
{
    using namespace std;
    using namespace vmc;

    MPIVMC::Init(); // make this usable with a MPI-compiled library

    // --- Create objects to move into VMC (don't use afterwards!)

    // Some trial wave functions
    auto psi1 = make_unique<Gaussian1D1POrbital>(1.2);
    auto psi2 = make_unique<Gaussian1D1POrbital>(0.5);
    auto psi3 = make_unique<Gaussian1D1POrbital>(1.0);
    auto psi4 = make_unique<QuadrExponential1D1POrbital>(-0.5, 1.0);

    // Hamiltonian for each wave function (keep in mind that the kinetic energy is strictly bound to it)
    // We use the harmonic oscillator with w=1
    auto ham1 = make_unique<HarmonicOscillator1D1P>(1.);
    auto ham2 = make_unique<HarmonicOscillator1D1P>(1.);
    auto ham3 = make_unique<HarmonicOscillator1D1P>(1.);
    auto ham4 = make_unique<HarmonicOscillator1D1P>(1.);


    // --- Create the VMC objects

    VMC vmc1(move(psi1), move(ham1));
    VMC vmc2(move(psi2), move(ham2));
    VMC vmc3(move(psi3), move(ham3));
    VMC vmc4(move(psi4), move(ham4));


    // --- Evaluate the energies

    cout << endl << " - - - EVALUATION OF ENERGY - - - " << endl << endl;

    const int E_NMC = 100000l; // MC samplings to use for computing the energy
    double energy[4]; // energy
    double d_energy[4]; // energy error bar

    // Case 1
    cout << "-> psi1: " << endl;
    vmc1.computeEnergy(E_NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    // Case 2
    cout << "-> psi2: " << endl;
    vmc2.computeEnergy(E_NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    // Case 3
    cout << "-> psi3: " << endl;
    vmc3.computeEnergy(E_NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    // Case 4
    cout << "-> psi4: " << endl;
    vmc4.computeEnergy(E_NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    MPIVMC::Finalize();

    return 0;
}
