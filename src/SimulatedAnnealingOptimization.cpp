#include "vmc/SimulatedAnnealingOptimization.hpp"

#include "vmc/MPIVMC.hpp"

namespace vmc
{
namespace vmc_siman
{
    WaveFunction * wf;
    Hamiltonian * H;
    int Nmc;
    mci::MCI * mci;

    // Simulated Annealing target function parameters
    double iota;   // factor that applies to the energy, for the target function
    double kappa;   // factor that applies to the energy's standard deviation, for the target function
    double lambda;  // normalization factor, for the target function

    // simulated annealing parameters
    gsl_siman_params_t params; // = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};


    inline double simanTarget(void *xp){
        using namespace std;

        // extract the variational parameters and apply them to the wf
        const auto * const x = (static_cast<double *>(xp));
        // apply the parameters to the wf
        wf->setVP(x);

        // compute the energy and its standard deviation
        double energy[4]; // energy
        double d_energy[4]; // energy error bar
        MPIVMC::Integrate(mci, Nmc, energy, d_energy, true, true);

        // compute the normalization
        const double norm = sqrt(std::inner_product(x, x+wf->getNVP(), x, 0.))/wf->getNVP();

        // return the target function
        return iota * energy[0] + kappa * d_energy[0] + lambda * norm;
    }

    inline void simanStep(const gsl_rng * r, void * xp, double step_size){
        // extract old variational parameters
        auto * const x = (static_cast<double *>(xp));

        // how many parameters should be changed? a random number between 1 and 1 + 2*log(N)
        const double numVariablesToChange = 1. + 2. * gsl_rng_uniform(r) * log( wf->getNVP() );
        const double rdThreshold = numVariablesToChange / wf->getNVP();

        // compute new variational parameters
        const double double_step_size = 2. * step_size;
        bool flag_change = false; // we want at least one variable changed
        do {
            for (int i=0; i<wf->getNVP(); ++i){
                if (gsl_rng_uniform(r)<rdThreshold) {
                    x[i] = x[i] + (gsl_rng_uniform(r)-0.5) * double_step_size;
                    flag_change = true;
                }
            }
        } while(!flag_change);
    }

    inline double simanDistance(void *xp, void *yp){
        // extract the two variatonal parameters of which I must compute the distance
        const auto * const x = (static_cast<double *>(xp));
        const auto * const y = (static_cast<double *>(yp));

        // compute the distance
        double dist = 0.;
        for (int i=0; i<wf->getNVP(); ++i){
            dist += (y[i] - x[i]) * (y[i] - x[i]);
        }
        return sqrt(dist);
    }

    inline void simanPrint(void *xp) {
        const auto * const vp = static_cast<double *>(xp);
        std::cout << "  [";
        for (int i=0; i<wf->getNVP(); ++i) {
            std::cout << " " << vp[i] << " ";
        }
        std::cout << "]  " << std::endl;
    }

} // namespace vmc_siman


SimulatedAnnealingOptimization::SimulatedAnnealingOptimization(WaveFunction * wf, Hamiltonian * H, const int &Nmc, mci::MCI * mci, const double &iota, const double &kappa, const double &lambda, gsl_siman_params_t &params)
    : WFOptimization(wf, H, mci)
{
    vmc_siman::wf = wf;
    vmc_siman::H = H;
    vmc_siman::Nmc = Nmc;
    vmc_siman::mci = mci;
    vmc_siman::iota = iota;
    vmc_siman::kappa = kappa;
    vmc_siman::lambda = lambda;
    vmc_siman::params = params;
}

void SimulatedAnnealingOptimization::optimizeWF()
{
    // random generator used by the simulated annealing
    gsl_rng_env_setup();
    const auto * T = gsl_rng_default;
    auto * r = gsl_rng_alloc(T);

    // set the initial parameters
    double vp[_wf->getNVP()];
    _wf->getVP(vp);

    // run the simulated annealing algorithm
    gsl_siman_solve(r, vp, vmc_siman::simanTarget, vmc_siman::simanStep, vmc_siman::simanDistance, vmc_siman::simanPrint, nullptr, nullptr, nullptr, _wf->getNVP()*sizeof(double), vmc_siman::params);

    // set the optimal variational parameters
    _wf->setVP(vp);

    // free resources
    gsl_rng_free(r);
}
} // namespace vmc