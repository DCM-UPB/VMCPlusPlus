#include "vmc/NMSimplexOptimization.hpp"
#include "vmc/MPIVMC.hpp"

#include <gsl/gsl_multimin.h>

namespace vmc
{

struct vmc_nms
{
    WaveFunction * wf;
    Hamiltonian * H;
    mci::MCI * mci;
    int Nmc;
    double iota;
    double kappa;
    double lambda;
    double rstart;
    double rend;
    size_t max_n_iter;

    void initFromOptimizer(NMSimplexOptimization * wfopt)
    {
        wf = wfopt->getWF();
        H = wfopt->getH();
        mci = wfopt->getMCI();
        Nmc = wfopt->getNmc();
        iota = wfopt->getIota();
        kappa = wfopt->getKappa();
        lambda = wfopt->getLambda();
        rstart = wfopt->getRStart();
        rend = wfopt->getREnd();
        max_n_iter = wfopt->getMaxNIter();
    }
};

double vmc_cost(const gsl_vector * v, void * params)
{
    WaveFunction * const wf = (static_cast<struct vmc_nms *>(params))->wf;
    Hamiltonian * const H = (static_cast<struct vmc_nms *>(params))->H;
    mci::MCI * const mci = (static_cast<struct vmc_nms *>(params))->mci;
    const int Nmc = (static_cast<struct vmc_nms *>(params))->Nmc;
    const double iota = (static_cast<struct vmc_nms *>(params))->iota;
    const double kappa = (static_cast<struct vmc_nms *>(params))->kappa;
    const double lambda = (static_cast<struct vmc_nms *>(params))->lambda;

    const auto nvp = static_cast<size_t>(wf->getNVP());
    double vpar[nvp];
    // apply the parameters to the wf
    for (size_t i = 0; i < nvp; ++i) {
        vpar[i] = gsl_vector_get(v, i);
    }
    wf->setVP(vpar);

    // compute the energy and its standard deviation
    double energy[4]; // energy
    double d_energy[4]; // energy error bar
    MPIVMC::Integrate(mci, Nmc, energy, d_energy, true, true);

    // compute the normalization
    double norm = 0.;
    for (size_t i = 0; i < nvp; ++i) {
        const double vi = gsl_vector_get(v, i);
        norm += vi*vi;
    }
    norm = sqrt(norm)/nvp;

    // return the cost function
    return iota*energy[0] + kappa*d_energy[0] + lambda*norm;
};


void NMSimplexOptimization::optimizeWF()
{
    const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer * s = nullptr;
    gsl_vector * ss, * x;
    gsl_multimin_function minex_func;

    int myrank = MPIVMC::MyRank();

    vmc_nms w{};
    w.initFromOptimizer(this);

    // Starting point
    const auto nvp = static_cast<size_t>(_wf->getNVP());
    double vpar[nvp];
    _wf->getVP(vpar);

    x = gsl_vector_alloc(nvp);
    for (size_t i = 0; i < nvp; ++i) {
        gsl_vector_set(x, i, vpar[i]);
    }

    // Set initial step sizes to 1
    ss = gsl_vector_alloc(nvp);
    gsl_vector_set_all(ss, _rstart);

    // Initialize method and iterate
    minex_func.n = nvp;
    minex_func.f = vmc_cost;
    minex_func.params = &w;

    s = gsl_multimin_fminimizer_alloc(T, nvp);
    gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

    size_t iter = 0;
    int status;
    do {
        status = gsl_multimin_fminimizer_iterate(s);

        if (status != 0) { break; }

        double size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, _rend);

        if (myrank == 0) {
            if (status == GSL_SUCCESS) {
                std::cout << "converged to minimum at" << std::endl;
            }
            std::cout << iter << " f() = " << s->fval << " size = " << size << std::endl;
        }
        ++iter;
    } while (status == GSL_CONTINUE && (_max_n_iter <= 0 || iter < _max_n_iter));

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
}
} // namespace vmc