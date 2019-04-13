#include <cmath>
#include <iostream>

#include "nfm/LogManager.hpp"
#include "vmc/Hamiltonian.hpp"
#include "vmc/MPIVMC.hpp"
#include "vmc/VMC.hpp"

#include "../common/ExampleFunctions.hpp"

int main()
{
    using namespace std;
    using namespace vmc;

    MPIVMC::Init(); // make this usable with a MPI-compiled library

    // Declare some trial wave functions
    auto * psi = new QuadrExponential1D1POrbital(-0.5, 1.0, true); // enable variational deriv

    // Declare an Hamiltonian
    // We use the harmonic oscillator with w=1 and w=2
    const double w = 1.;
    auto * ham = new HarmonicOscillator1D1P(w, psi);

    nfm::LogManager::setLoggingOn(false); // use this to enable log printout (true would mean verbose mode)

    cout << endl << " - - - WAVE FUNCTION OPTIMIZATION - - - " << endl << endl;

    const int NMC = 10000l; // MC samplings to use for computing the energy
    double energy[4]; // energy
    double d_energy[4]; // energy error bar
    double vp[psi->getNVP()];


    auto * vmc = new VMC(psi, ham);
    cout << "-> ham:    w = " << w << endl << endl;

    cout << "   Initial Wave Function parameters:" << endl;
    psi->getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Starting energy:" << endl;
    vmc->computeEnergy(NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    // settings for better performance
    vmc->getMCI()->setNfindMRT2Iterations(10);
    vmc->getMCI()->setNdecorrelationSteps(1000);

    cout << "   Optimization . . ." << endl;
    vmc->stochasticReconfigurationOptimization(NMC, 0.01, true);
    cout << "   . . . Done!" << endl << endl;

    cout << "   Optimized Wave Function parameters:" << endl;
    psi->getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Optimized energy:" << endl;
    vmc->computeEnergy(NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl << endl;


    delete vmc;
    delete ham;
    delete psi;

    MPIVMC::Finalize();

    return 0;
}
