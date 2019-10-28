#include <cmath>
#include <cassert>
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

    const bool verbose = false;
    const double w = 1.0; // harmonic oscillator strength

    // uncomment to check a range
    //const std::vector<double> ps{0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
    //                             0.925, 0.95, 0.96, 0.97, 0.98, 0.99, 1.01, 1.02, 1.03, 1.04, 1.05, 1.075,
    //                             1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};

    // uncomment to check only one position
    const std::vector<double> ps{1.2};

    std::vector<double> ens_vmc;
    std::vector<double> ens_ana;
    std::vector<double> diff_ens;
    std::vector<double> reldiff_grad;
    std::vector<double> relerr_grad;

    for (const double p : ps) { // gaussian wf variational parameter

        // --- Create the VMC object
        const int nskip = 2; // compute energy/grad only every second MC step
        const int blocksize = 8; // we use fixed block size for performance and memory
        VMC vmc(make_unique<ConstNormGaussian1D1POrbital>(p), make_unique<HarmonicOscillator1D1P>(w), nskip, blocksize);
        vmc.getMCI().setTrialMove(mci::SRRDType::Gaussian);
        vmc.getMCI().setSeed(1337 + 42*myrank); // we need to use a fixed seed such to make sure that the noisy asserts will always pass

        // --- Evaluate the energies
        const int E_NMC = 1024*1024; // MC samplings to use for computing the energy (kept low to allow running inside valgrind)

        // compute energy and gradient via VMC sampling
        EnergyGradientTargetFunction gradfun(vmc, E_NMC, E_NMC, true, 0.); // create gradient target function
        std::vector<double> x0(1);
        vmc.getVP(x0.data()); // set x0 from wf VP
        nfm::NoisyGradient grad(1);
        const auto en = gradfun.fgrad(x0, grad);

        // For the p-parametrized gaussian
        // Psi(x) ~ sqrt(p) e(-0.5*p^2*x^2)
        // the energy should be E(p) = (w^2 + p^4)/(4 p^2)
        // and the gradient d/dp E(p) = (-w^2 + p^4)/(2 p^3)
        const double en_ana = (w*w + p*p*p*p)/(4.*p*p);
        const double grad_ana = (w*w - p*p*p*p)/(2.*p*p*p);

        const double en_diff = fabs(en.val - en_ana);
        const double grad_diff = fabs(grad.val[0] - grad_ana);

        // report, if enabled
        if (myrank == 0 && verbose) {
            cout << endl;
            cout << "Energy (VMC)     = " << en.val << " +- " << en.err << " (rel: " << en.err/fabs(en_ana) << ")" << endl;
            cout << "Energy (ANA)     = " << en_ana << endl;
            cout << "Energy Diff      = " << en_diff << endl;
            cout << "Gradient (VMC)   = " << grad.val[0] << " +- " << grad.err[0] << " (rel: " << grad.err[0]/fabs(grad_ana) << ")" << endl;
            cout << "Gradient (ANA) = " << grad_ana << endl;
            cout << "Gradient Diff  = " << grad_diff << endl;
            cout << endl;
        }

        // use 3-sigma asserts
        assert(en_diff < 3.*en.err);
        assert(grad_diff < 3.*grad.err[0]);

        ens_vmc.push_back(en.val);
        ens_ana.push_back(en_ana);
        diff_ens.push_back(fabs(en.val-en_ana));
        reldiff_grad.push_back(fabs(grad_diff/grad_ana));
        relerr_grad.push_back(fabs(grad.err[0]/grad_ana));
    }

    if (myrank == 0 && verbose) { // will report results over range
        cout << "Curve results: " << endl;
        cout << "p   E        E_ana    E - E_ana  |(g-g_ana)/g_ana| |dg/g_ana|" << endl;
        for (size_t i = 0; i < ps.size(); ++i) {
            cout << ps[i] << " " << ens_vmc[i] << " " << ens_ana[i] << " ";
            cout << diff_ens[i] << " " << reldiff_grad[i] << " " << relerr_grad[i] << endl;
        }
    }

    MPIVMC::Finalize();

    return 0;
}
