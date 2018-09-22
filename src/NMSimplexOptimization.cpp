#include "NMSimplexOptimization.hpp"
#include "MPIVMC.hpp"

struct vmc_workspace
{
    WaveFunction * wf;
    Hamiltonian * H;
    MCI * mci;
    long Nmc;
    double iota;
    double kappa;
    double lambda;

    void initFromOptimizer(NMSimplexOptimization * wfopt)
    {
        wf = wfopt->getWF();
        H = wfopt->getH();
        mci = wfopt->getMCI();
        Nmc = wfopt->getNmc();
        iota = wfopt->getIota();
        kappa = wfopt->getKappa();
        lambda = wfopt->getLambda();
    }
};

double vmc_cost(const gsl_vector *v, void *params)
{
    WaveFunction * const wf = ((struct vmc_workspace *)params)->wf;
    Hamiltonian * const H = ((struct vmc_workspace *)params)->H;
    MCI * const mci = ((struct vmc_workspace *)params)->mci;
    const long Nmc = ((struct vmc_workspace *)params)->Nmc;
    const double iota = ((struct vmc_workspace *)params)->iota;
    const double kappa = ((struct vmc_workspace *)params)->kappa;
    const double lambda = ((struct vmc_workspace *)params)->lambda;

    double vpar[wf->getNVP()];
    // apply the parameters to the wf
    for (int i=0; i<wf->getNVP(); ++i) {
        vpar[i] = gsl_vector_get(v, i);
    }
    wf->setVP(&vpar[0]);

    // compute the energy and its standard deviation
    double energy[4]; // energy
    double d_energy[4]; // energy error bar
    mci->clearSamplingFunctions(); mci->addSamplingFunction(wf);
    mci->clearObservables(); mci->addObservable(H);
    MPIVMC::Integrate(mci, Nmc, energy, d_energy, true, true);

    // compute the normalization
    double norm = 0.;
    for (int i=0; i<wf->getNVP(); ++i)
        norm += pow(gsl_vector_get(v, i), 2);
    norm = sqrt(norm)/wf->getNVP();

    // return the cost function
    return iota * energy[0] + kappa * d_energy[0] + lambda * norm;
};


void NMSimplexOptimization::optimizeWF()
{

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  vmc_workspace w;
  w.initFromOptimizer(this);

  // Starting point
  int npar = _wf->getNVP();
  double vpar[npar];
  _wf->getVP(&vpar[0]);

  x = gsl_vector_alloc(npar);
  for (int i=0; i<npar; ++i) {
      gsl_vector_set(x, i, vpar[i]);
  }

  // Set initial step sizes to 1
  ss = gsl_vector_alloc (npar);
  gsl_vector_set_all(ss, 1.0);

  // Initialize method and iterate
  minex_func.n = npar;
  minex_func.f = vmc_cost;
  minex_func.params = &w;

  s = gsl_multimin_fminimizer_alloc(T, npar);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status)
        break;

      size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, 0.01);

      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }

      printf ("%5zu f() = %7.3f size = %.3f\n", iter, s->fval, size);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}
