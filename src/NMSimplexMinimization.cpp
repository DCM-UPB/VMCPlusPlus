#include "vmc/NMSimplexMinimization.hpp"
#include "vmc/MPIVMC.hpp"

#include <gsl/gsl_multimin.h>

namespace vmc
{

struct vmc_nms
{
    VMC &vmc;
    int Nmc;
    double iota;
    double kappa;
    double lambda;
    double rstart;
    double rend;
    size_t max_n_iter;

    vmc_nms(VMC &vmc, NMSimplexMinimization &nms): vmc(vmc)
    {
        Nmc = nms.getNmc();
        iota = nms.getIota();
        kappa = nms.getKappa();
        lambda = nms.getLambda();
        rstart = nms.getRStart();
        rend = nms.getREnd();
        max_n_iter = nms.getMaxNIter();
    }
};

double vmc_cost(const gsl_vector * v, void * params)
{
    VMC &vmc = (static_cast<struct vmc_nms *>(params))->vmc;
    const int Nmc = (static_cast<struct vmc_nms *>(params))->Nmc;
    const double iota = (static_cast<struct vmc_nms *>(params))->iota;
    const double kappa = (static_cast<struct vmc_nms *>(params))->kappa;
    const double lambda = (static_cast<struct vmc_nms *>(params))->lambda;

    const auto nvp = static_cast<size_t>(vmc.getNVP());
    double vpar[nvp];
    // apply the parameters to the wf
    for (size_t i = 0; i < nvp; ++i) {
        vpar[i] = gsl_vector_get(v, i);
    }
    vmc.setVP(vpar);

    // compute the energy and its standard deviation
    double energy[4]; // energy
    double d_energy[4]; // energy error bar
    vmc.computeEnergy(Nmc, energy, d_energy, true, true);

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


void NMSimplexMinimization::minimizeEnergy(VMC &vmc)
{
    const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer * s = nullptr;
    gsl_vector * ss, * x;
    gsl_multimin_function minex_func;

    int myrank = MPIVMC::MyRank();

    vmc_nms w(vmc, *this);

    // Starting point
    const auto nvp = static_cast<size_t>(vmc.getNVP());
    double vpar[nvp];
    vmc.getVP(vpar);

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