#include <cmath>
#include <iostream>
#include <memory>

#include "vmc/Hamiltonian.hpp"
#include "vmc/MPIVMC.hpp"
#include "vmc/VMC.hpp"

#include "TestVMCFunctions.hpp"


int main()
{
    using namespace std;
    using namespace vmc;

    MPIVMC::Init();

    const int NMC = 4000l;

    auto gauss = std::make_unique<Gaussian1D1POrbital>(1.2); // will be moved from
    auto harm_osc = std::make_unique<HarmonicOscillator1D1P>(1., gauss.get()); // will be moved from
    VMC vmc(std::move(gauss), std::move(harm_osc));

    cout << endl << " - - - EVALUATION OF ENERGY - - - " << endl << endl;
    double b1;
    vmc.getWF().getVP(&b1);
    cout << "Wave Function b     = " << b1 << endl;
    vmc.getMCI().setIRange(-10., 10.);
    double energy[4];
    double d_energy[4];
    vmc.computeEnergy(NMC, energy, d_energy, false, false);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    cout << endl << " - - - ONE-DIMENSIONAL MINIMIZATION - - - " << endl << endl;
    double b;
    vmc.getWF().getVP(&b);
    cout << "Wave Function b     = " << b << endl;
    cout << "Conjugate Gradient Minimization ... " << endl;
    vmc.conjugateGradientOptimization(NMC, 4*NMC);
    vmc.getWF().getVP(&b);
    cout << "Wave Function b     = " << b << endl << endl;
    vmc.computeEnergy(NMC, energy, d_energy);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    MPIVMC::Finalize();

    return 0;
}
