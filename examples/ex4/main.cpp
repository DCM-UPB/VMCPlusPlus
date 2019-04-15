#include <cmath>
#include <gsl/gsl_siman.h>
#include <iostream>

#include "vmc/MPIVMC.hpp"
#include "vmc/VMC.hpp"
#include "vmc/SimulatedAnnealingMinimization.hpp"

#include "../common/ExampleFunctions.hpp"

int main()
{
    using namespace std;
    using namespace vmc;

    MPIVMC::Init(); // make this usable with a MPI-compiled library

    // Setup VMC
    auto psi = make_unique<QuadrExponential1D1POrbital>(-0.5, 1.0, true); // enable variational deriv
    const double w = 1.; // We use the harmonic oscillator with w=1
    auto ham = make_unique<HarmonicOscillator1D1P>(1.);
    VMC vmc(move(psi), move(ham));


    cout << endl << " - - - WAVE FUNCTION OPTIMIZATION - - - " << endl << endl;

    const int NMC = 10000l; // MC samplings to use for computing the energy
    double energy[4]; // energy
    double d_energy[4]; // energy error bar
    double vp[vmc.getNVP()];


    cout << "-> ham:    w = " << w << endl << endl;

    cout << "   Initial Wave Function parameters:" << endl;
    vmc.getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Starting energy:" << endl;
    vmc.computeEnergy(NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    cout << "   Optimization . . ." << endl;

    // settings for performance
    vmc.getMCI().setNfindMRT2Iterations(10);
    vmc.getMCI().setNdecorrelationSteps(1000);

    // simulated annealing parameters
    int N_TRIES = 20;
    int ITERS_FIXED_T = 20;
    double STEP_SIZE = 0.1;
    double K = 0.1;
    double T_INITIAL = 1.0;
    double MU_T = 1.3;
    double T_MIN = 0.01;
    gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};

    SimulatedAnnealingMinimization siman(NMC, 1., 0.1, 0., params);
    siman.minimizeEnergy(vmc);

    cout << "   . . . Done!" << endl << endl;

    cout << "   Optimized Wave Function parameters:" << endl;
    vmc.getVP(vp);
    cout << "       a = " << vp[0] << endl;
    cout << "       b = " << vp[1] << endl;

    cout << "   Optimized energy:" << endl;
    vmc.computeEnergy(NMC, energy, d_energy);
    cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl << endl;


    MPIVMC::Finalize();

    return 0;
}
