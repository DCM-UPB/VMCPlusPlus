#include <cmath>
#include <iostream>

#include "nfm/Adam.hpp"
#include "nfm/LogManager.hpp"
#include "vmc/MPIVMC.hpp"
#include "vmc/VMC.hpp"
#include "vmc/EnergyMinimization.hpp"

#include "../common/ExampleFunctions.hpp"

int main()
{
    using namespace std;
    using namespace vmc;

    MPIVMC::Init(); // make this usable with a MPI-compiled library


    // --- Create objects to move into VMC (don't use afterwards!)

    // Create some trial wave functions
    auto psi1 = make_unique<QuadrExponential1D1POrbital>(-0.5, 1.0, true); // enable variational derivatives
    auto psi2 = make_unique<QuadrExponential1D1POrbital>(-0.5, 1.0, true); // we create an equivalent psi2

    // Create corresponding Hamiltonians
    // We use the harmonic oscillator with w=1 and w=2
    const double w1 = 1.;
    auto ham1 = make_unique<HarmonicOscillator1D1P>(w1);
    const double w2 = 2.;
    auto ham2 = make_unique<HarmonicOscillator1D1P>(w2);

    // --- Create the VMC objects

    VMC vmc1(move(psi1), move(ham1));
    VMC vmc2(move(psi2), move(ham2));


    if (MPIVMC::MyRank() == 0) {
        nfm::LogManager::setLoggingOn(false); // use this to enable log printout (true would mean verbose mode)
    }

    cout << endl << " - - - WAVE FUNCTION OPTIMIZATION - - - " << endl << endl;

    const int E_NMC = 4000l; // MC samplings to use for computing the energy
    const int G_NMC = 10000l; // MC samplings to use for computing the energy & gradient
    double energy[4]; // energy
    double d_energy[4]; // energy error bar
    double vp[vmc1.getNVP()];


    // Case 1
    cout << "-> ham1:    w = " << w1 << endl << endl;

    const double stepSize = 0.1; // in this case we should set a larger ADAM step size than default (0.001)

    cout << "   Initial Wave Function parameters:" << endl;
    vmc1.getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Starting energy:" << endl;
    vmc1.computeEnergy(E_NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    cout << "   Optimization . . ." << endl;

    // settings for better performance
    vmc1.getMCI().setNfindMRT2Iterations(10);
    vmc1.getMCI().setNdecorrelationSteps(1000);

    nfm::Adam adam(vmc1.getNVP(), true /* use averaging for final parameters */, stepSize);
    adam.setMaxNConstValues(10); // stop after 10 constant values (within error)
    minimizeEnergy(vmc1, adam, E_NMC, G_NMC, false /* don't use/calculate gradient error */, 0.01 /* parameter regularization */);
    cout << "   . . . Done!" << endl << endl;

    cout << "   Optimized Wave Function parameters:" << endl;
    vmc1.getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Optimized energy:" << endl;
    vmc1.computeEnergy(E_NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl << endl;



    // Case 2
    cout << "-> ham2:    w = " << w2 << endl << endl;

    cout << "   Initial Wave Function parameters:" << endl;
    vmc2.getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Starting energy:" << endl;
    vmc2.computeEnergy(E_NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    cout << "   Optimization . . ." << endl;

    // settings for better performance
    vmc2.getMCI().setNfindMRT2Iterations(10);
    vmc2.getMCI().setNdecorrelationSteps(1000);

    minimizeEnergy(vmc2, adam, E_NMC, G_NMC, false, 0.01);
    cout << "   . . . Done!" << endl << endl;

    cout << "   Optimized Wave Function parameters:" << endl;
    vmc2.getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Optimized energy:" << endl;
    vmc2.computeEnergy(E_NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;


    MPIVMC::Finalize();

    return 0;
}
