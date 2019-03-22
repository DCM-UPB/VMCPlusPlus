#include "vmc/NMSimplexOptimization.hpp"
#include "vmc/MPIVMC.hpp"

#include <gsl/gsl_multimin.h>

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
    int max_n_iter;

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

double vmc_cost(const gsl_vector *v, void *params)
{
    WaveFunction * const wf = (static_cast<struct vmc_nms *>(params))->wf;
    Hamiltonian * const H = (static_cast<struct vmc_nms *>(params))->H;
    mci::MCI * const mci = (static_cast<struct vmc_nms *>(params))->mci;
    const int Nmc = (static_cast<struct vmc_nms *>(params))->Nmc;
    const double iota = (static_cast<struct vmc_nms *>(params))->iota;
    const double kappa = (static_cast<struct vmc_nms *>(params))->kappa;
    const double lambda = (static_cast<struct vmc_nms *>(params))->lambda;

    double vpar[wf->getNVP()];
    // apply the parameters to the wf
    for (int i=0; i<wf->getNVP(); ++i) {
        vpar[i] = gsl_vector_get(v, i);
    }
    wf->setVP(vpar);

    // compute the energy and its standard deviation
    double energy[4]; // energy
    double d_energy[4]; // energy error bar
    mci->clearSamplingFunctions(); mci->addSamplingFunction(*wf);
    mci->clearObservables(); mci->addObservable(*H);
    MPIVMC::Integrate(mci, Nmc, energy, d_energy, true, true);

    // compute the normalization
    double norm = 0.;
    for (int i=0; i<wf->getNVP(); ++i) {
        const double vi = gsl_vector_get(v, i);
        norm += vi*vi;
    }
    norm = sqrt(norm)/wf->getNVP();

    // return the cost function
    return iota * energy[0] + kappa * d_energy[0] + lambda * norm;
};


void NMSimplexOptimization::optimizeWF()
{
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = nullptr;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    int myrank = MPIVMC::MyRank();

    vmc_nms w{};
    w.initFromOptimizer(this);

    // Starting point
    int npar = _wf->getNVP();
    double vpar[npar];
    _wf->getVP(vpar);

    x = gsl_vector_alloc(npar);
    for (int i=0; i<npar; ++i) {
        gsl_vector_set(x, i, vpar[i]);
    }

    // Set initial step sizes to 1
    ss = gsl_vector_alloc (npar);
    gsl_vector_set_all(ss, _rstart);

    // Initialize method and iterate
    minex_func.n = npar;
    minex_func.f = vmc_cost;
    minex_func.params = &w;

    s = gsl_multimin_fminimizer_alloc(T, npar);
    gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

    int iter = 0;
    int status;
    do {
        status = gsl_multimin_fminimizer_iterate(s);

        if (status != 0) { break; }

        double size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, _rend);

        if (myrank==0) {
            if (status == GSL_SUCCESS) {
                std::cout << "converged to minimum at" << std::endl;
            }
            std::cout << iter << " f() = " << s->fval << " size = " << size << std::endl;
        }
        ++iter;
    } while (status == GSL_CONTINUE && (iter < _max_n_iter));

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
}
