#include <cmath>
#include <iostream>
#include <memory>

#include "vmc/MPIVMC.hpp"
#include "vmc/VMC.hpp"
#include "vmc/EnergyGradientTargetFunction.hpp"

#include "../common/TestVMCFunctions.hpp" // used WFs / Hamiltonian


int main()
{
    using namespace std;
    using namespace vmc;

    const int myrank = MPIVMC::Init(); // make this usable with a MPI-compiled library

    // uncomment to check a range
    //const std::vector<double> p0s{0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
    //                              0.925, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.075,
    //                              1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};

    // uncomment to check only one position
    const std::vector<double> p0s{1.01};

    std::vector<double> reldiffs_num;
    std::vector<double> enls;

    for (const double p : p0s) {

        // --- Create the VMC objects

        const double w = 1.0; // harmonic oscillator strength
        const double p0 = p; // gaussian wf variational parameter
        const double dp = 0.005; // change for numerical gradient calculation
        const int nskip = 2; // compute energy/grad only every second MC step
        const int blocksize = 8; // we use fixed block size for performance and memory
        VMC vmc_l(make_unique<ConstNormGaussian1D1POrbital>(p0), make_unique<HarmonicOscillator1D1P>(w), nskip, blocksize);
        VMC vmc_m(make_unique<ConstNormGaussian1D1POrbital>(p0 + dp), make_unique<HarmonicOscillator1D1P>(w), nskip, blocksize);
        VMC vmc_r(make_unique<ConstNormGaussian1D1POrbital>(p0 + dp + dp), make_unique<HarmonicOscillator1D1P>(w), nskip, blocksize);
        vmc_l.getMCI().setTrialMove(mci::SRRDType::Gaussian);
        vmc_m.getMCI().setTrialMove(mci::SRRDType::Gaussian);
        vmc_r.getMCI().setTrialMove(mci::SRRDType::Gaussian);

        // --- Evaluate the energies

        const int E_NMC = 1024*1024*512; // MC samplings to use for computing the energy

        // first compute energy and gradient at position 1
        EnergyGradientTargetFunction gradfun(vmc_m, E_NMC, E_NMC, true, 0.); // create gradient target function
        std::vector<double> x0(1);
        vmc_m.getVP(x0.data()); // set x0 from wf VP
        nfm::NoisyGradient grad(1);
        const auto en = gradfun.fgrad(x0, grad);

        // Now the energy at left position
        double ens_l[4]; // energies
        double d_ens_l[4]; // energies error bar
        vmc_l.computeEnergy(E_NMC, ens_l, d_ens_l);

        // Now the energy at right position
        double ens_r[4]; // energies
        double d_ens_r[4]; // energies error bar
        vmc_r.computeEnergy(E_NMC, ens_r, d_ens_r);

        // And the numeric gradient (negative direction)
        nfm::NoisyGradient grad_num(1);
        grad_num.val[0] = -(ens_r[0] - ens_l[0])/(dp + dp);
        grad_num.err[0] = sqrt(pow(d_ens_r[0], 2) + pow(d_ens_l[0], 2))/(dp + dp);

        if (myrank == 0) {
            cout << endl;
            cout << "Energy Left = " << ens_l[0] << " +- " << d_ens_l[0] << " (rel: " << d_ens_l[0]/fabs(ens_l[0]) << ")" << endl;
            cout << "Energy Mid  = " << en.val << " +- " << en.err << " (rel: " << en.err/fabs(en.val) << ")" << endl;
            cout << "Energy Right = " << ens_r[0] << " +- " << d_ens_r[0] << " (rel: " << d_ens_r[0]/fabs(ens_r[0]) << ")" << endl;
            cout << "Gradient (ana) = " << grad.val[0] << " +- " << grad.err[0] << " (rel: " << grad.err[0]/fabs(grad.val[0]) << ")" << endl;
            cout << "Gradient (num) = " << grad_num.val[0] << " +- " << grad_num.err[0] << " (rel: " << grad_num.err[0]/fabs(grad_num.val[0]) << ")" << endl;
            cout << endl;
        }

        enls.push_back(ens_l[0]);
        reldiffs_num.push_back(grad_num.err[0]/fabs(grad_num.val[0]));
    }

    if (myrank == 0 && p0s.size() > 1) { // will report if range was used
        cout << "Curve results: " << endl;
        cout << "p0    E    dg/abs(g)" << endl;
        for (size_t i = 0; i < p0s.size(); ++i) {
            cout << p0s[i] << " " << enls[i] << " " << reldiffs_num[i] << endl;
        }
    }

    // TODO:
    // Test should better compare MC-sampled to true analytical gradients,
    // instead of relying on biased and noisy finite difference gradients.
    // For w=1 harmonic oscillator and a normalized gaussian wave function
    // Psi(x) ~ sqrt(p) e(-0.5*p^2*x^2) , i.e. <Psi(p)|Psi(p)> is constant for all p
    // the energy should be E(p) = (1 + p^4)/(4 p^2)
    // and the gradient d/dp E(p) = (-1 + b^4)/(2 b^3)

    MPIVMC::Finalize();

    return 0;
}
