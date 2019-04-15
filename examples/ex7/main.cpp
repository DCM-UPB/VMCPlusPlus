#include <cmath>
#include <iostream>

#include "nfm/ConjGrad.hpp"
#include "vmc/Hamiltonian.hpp"
#include "vmc/VMC.hpp"

#include "vmc/MPIVMC.hpp" // this example requires MPI!

#include "../common/ExampleFunctions.hpp"

int main()
{
    using namespace std;
    using namespace vmc;

    int myrank = MPIVMC::Init();
    cout << "Hello from rank " << myrank << endl;

    // Setup VMC
    auto psi = make_unique<Gaussian1D1POrbital>(0.6); // we use a simple gaussian, but not with ground state parametrization
    const double w = 1.; // We use the harmonic oscillator with w=1
    auto ham = make_unique<HarmonicOscillator1D1P>(1.);
    VMC vmc(move(psi), move(ham));


    const int E_NMC = 100000l; // MC samplings to use for computing the energy
    double energy[4], energy_h[4]; // energy
    double d_energy[4], d_energy_h[4]; // energy error bar
    for (int i = 0; i < 4; ++i) {
        energy[i] = 0.;
        d_energy[i] = 0.;
    }

    // example of file out with MPI (but we disable it for the test below (100K Steps!!))
    auto obsfile = "obsfile" + std::to_string(myrank);
    auto wlkfile = "wlkfile" + std::to_string(myrank);
    vmc.getMCI().storeObservablesOnFile(obsfile, 1); // would print observables on every step
    vmc.getMCI().storeWalkerPositionsOnFile(wlkfile, 1); // would print walker positions on every step
    vmc.getMCI().clearObservableFile();
    vmc.getMCI().clearWalkerFile();

    if (myrank == 0) {
        cout << endl << " - - - EVALUATION OF ENERGY - - - " << endl << endl;
    }

    const int neval = 5;
    // First compute the energy with auto-mode findMRT2step/initialDecorr/blocking
    if (myrank == 0) {
        cout << "Computing the energy " << neval << " times, with auto-mode findMRT2step/initialDecorr (inconsistent time per CPU)." << endl;
    }
    for (int i = 0; i < neval; ++i) {
        vmc.computeEnergy(E_NMC, energy_h, d_energy_h);
        if (myrank == 0) {
            for (int j = 0; j < 4; ++j) {
                energy[j] += energy_h[j];
                d_energy[j] += d_energy_h[j]*d_energy_h[j];
            }
            cout << "Total Energy        = " << energy_h[0] << " +- " << d_energy_h[0] << endl;
        }
    }
    if (myrank == 0) {
        cout << "On average:" << endl;
        cout << "Total Energy        = " << energy[0]/neval << " +- " << sqrt(d_energy[0])/neval << endl;
        cout << "Potential Energy    = " << energy[1]/neval << " +- " << sqrt(d_energy[1])/neval << endl;
        cout << "Kinetic (PB) Energy = " << energy[2]/neval << " +- " << sqrt(d_energy[2])/neval << endl;
        cout << "Kinetic (JF) Energy = " << energy[3]/neval << " +- " << sqrt(d_energy[3])/neval << endl << endl;
    }
    for (int i = 0; i < 4; ++i) {
        energy[i] = 0.;
        d_energy[i] = 0.;
    }

    // Now fix the number of steps for findMRT2step/initialDecorr
    // we use a generous total amount of 10000 equilibration steps
    vmc.getMCI().setNfindMRT2Iterations(50);
    vmc.getMCI().setNdecorrelationSteps(5000);

    if (myrank == 0) {
        cout << "Computing the energy " << neval << " times, with fixed-mode findMRT2step/initialDecorr (consistent time per CPU)." << endl;
    }
    for (int i = 0; i < neval; ++i) {
        vmc.computeEnergy(E_NMC, energy_h, d_energy_h);
        if (myrank == 0) {
            for (int j = 0; j < 4; ++j) {
                energy[j] += energy_h[j];
                d_energy[j] += d_energy_h[j]*d_energy_h[j];
            }
            cout << "Total Energy        = " << energy_h[0] << " +- " << d_energy_h[0] << endl;
        }
    }
    if (myrank == 0) {
        cout << "On average:" << endl;
        cout << "Total Energy        = " << energy[0]/neval << " +- " << sqrt(d_energy[0])/neval << endl;
        cout << "Potential Energy    = " << energy[1]/neval << " +- " << sqrt(d_energy[1])/neval << endl;
        cout << "Kinetic (PB) Energy = " << energy[2]/neval << " +- " << sqrt(d_energy[2])/neval << endl;
        cout << "Kinetic (JF) Energy = " << energy[3]/neval << " +- " << sqrt(d_energy[3])/neval << endl << endl;
    }

    MPIVMC::Finalize();

    return 0;
}
