#include "StochasticReconfigurationOptimization.hpp"
#include "StochasticReconfigurationMCObservable.hpp"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>

namespace sropt_details {
    void vmc_workspace::initFromOptimizer(StochasticReconfigurationOptimization * wfopt)
    {
        wf = wfopt->getWF();
        H = wfopt->getH();
        mci = wfopt->getMCI();
        Nmc = wfopt->getNmc();
    }

    void vmc_integrate(vmc_workspace &w, const double * const vp, double * const obs, double * const dobs, const bool flag_grad = false)
    {
        // set the variational parameters given as input
        w.wf->setVP(vp);

        // set up the MC integrator
        w.mci->clearSamplingFunctions();
        w.mci->addSamplingFunction(w.wf);

        w.mci->clearObservables();
        w.mci->addObservable(w.H);
        if(flag_grad) w.mci->addObservable(new StochasticReconfigurationMCObservable(w.wf, w.H));

        // perform the integral and store the values
        w.mci->integrate(w.Nmc, obs, dobs);
    }

    void vmc_calcobs(vmc_workspace &w, const double * const vp, double &f, double &df, double * const grad_E = NULL , double * const dgrad_E = NULL)
    {
        const int nvp = w.wf->getNVP();
        const bool flag_grad = (grad_E && dgrad_E);

        double obs[flag_grad ? 4 + 2*nvp + nvp*nvp : 4];
        double dobs[flag_grad ? 4 + 2*nvp + nvp*nvp : 4];

        vmc_integrate(w, vp, obs, dobs, flag_grad);

        f = obs[0];
        df = dobs[0];

        if (flag_grad) {
            // create pointers for ease of use and readability
            double * const H = obs;
            double * const dH = dobs;
            double * const Oi = obs+4;
            double * const dOi = dobs+4;
            double * const HOi = obs+4+nvp;
            double * const dHOi = dobs+4+nvp;
            double * const OiOj = obs+4+2*nvp;
            double * const dOiOj = dobs+4+2*nvp;


            // --- compute direction (or gradient) to follow
            gsl_matrix * sij = gsl_matrix_alloc(nvp, nvp);
            gsl_matrix * rdsij = gsl_matrix_alloc(nvp, nvp);   // relative error, i.e. error/value
            for (int i=0; i<nvp; ++i){
                for (int j=0; j<nvp; ++j){
                    gsl_matrix_set(sij, i, j, OiOj[i*nvp + j] - Oi[i] * Oi[j]);
                    gsl_matrix_set(rdsij, i, j,
                                   (dOiOj[i*nvp + j] + abs(Oi[i]*Oi[j])*( (dOi[i]/Oi[i]) + (dOi[j]/Oi[j]) ))
                                   /  gsl_matrix_get(sij, i, j) );
                }
            }
            gsl_vector * fi = gsl_vector_alloc(nvp);
            gsl_vector * rdfi = gsl_vector_alloc(nvp);   // relative error, i.e. error/value
            for (int i=0; i<nvp; ++i){
                gsl_vector_set(fi, i, H[0]*Oi[i] - HOi[i]);
                gsl_vector_set(rdfi, i,
                               (abs(H[0]*Oi[i])*( (dH[0]/H[0]) + (dOi[i]/Oi[i]) ) + dHOi[i])
                               / gsl_vector_get(fi, i)  );
            }
            // invert matrix using SVD
            const double SVD_MIN = 1.0e-9;
            // matrix and vectors needed for the SVD
            gsl_matrix * V = gsl_matrix_alloc(nvp, nvp);
            gsl_vector * S = gsl_vector_alloc(nvp);
            gsl_vector * work = gsl_vector_alloc(nvp);
            // run the Single Value Decomposition
            gsl_linalg_SV_decomp(sij, V, S, work);
            // assemble the inverse matrix
            gsl_matrix * Isij = gsl_matrix_alloc(nvp, nvp);
            for (int i=0; i<nvp; ++i){
                for (int j=0; j<nvp; ++j){
                    gsl_matrix_set(Isij, i, j, 0.);
                }
            }
            for (int i=0; i<nvp; ++i){
                if (gsl_vector_get(S, i) > SVD_MIN * gsl_vector_get(S, 0)){
                    gsl_matrix_set(Isij, i, i,   1./gsl_vector_get(S, i)  );
                } else {
                    gsl_matrix_set(Isij, i, i,   0.  );
                }
            }
            gsl_matrix * mm = gsl_matrix_alloc(nvp, nvp);
            gsl_matrix_transpose(V);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Isij, V, 0.0, mm);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, sij, mm, 0.0, Isij);
            // --- finally find the direction to follow
            double foo;
            for (int i=0; i<nvp; ++i){
                grad_E[i] = 0.;
                dgrad_E[i] = 0.;
                for (int k=0; k<nvp; ++k){
                    foo = gsl_vector_get(fi, k)*gsl_matrix_get(Isij, k, i);
                    grad_E[i] -= foo;  // the minus sign is there beacuse the NoisyFunMin library will follow 'minus the gradient'
                    dgrad_E[i] += abs(foo) * ( gsl_vector_get(rdfi, k) + gsl_matrix_get(rdsij, k, i) );  // not correct, just a rough estimation
                }
            }

            // free resources
            gsl_matrix_free(mm);
            gsl_vector_free(work);
            gsl_vector_free(S);
            gsl_matrix_free(V);
            gsl_matrix_free(Isij);
            gsl_vector_free(rdfi);
            gsl_vector_free(fi);
            gsl_matrix_free(rdsij);
            gsl_matrix_free(sij);
        }
    }

    // NoisyFunctionWithGradient implementation
    void fval(vmc_workspace &w, const double *vp, double &f, double &df){
        vmc_calcobs(w, vp, f, df);
    }

    void grad(vmc_workspace &w, const double *vp, double *grad_E, double *dgrad_E){
        double f, df;
        vmc_calcobs(w, vp, f, df, grad_E, dgrad_E);
    }

    void fgrad(vmc_workspace &w, const double *vp, double &f, double &df, double *grad_E, double *dgrad_E){
        vmc_calcobs(w, vp, f, df, grad_E, dgrad_E);
    }


    // for GSL version
    void vmc_vpar(const gsl_vector * const v, double * const vpar)
    {
        for (int i=0; i<(int)v->size; ++i) {
            vpar[i] = gsl_vector_get(v, i);
        }
    }

    void vmc_cost_grad(const gsl_vector *v, void *params, double *f, gsl_vector *df)
    {
        if (!f && !df) return; // both NULL

        WaveFunction * const wf = ((struct vmc_workspace *)params)->wf;

        double vpar[wf->getNVP()];
        vmc_vpar(v, vpar);

        double err;
        if (f && !df) {
            fval(*(struct vmc_workspace *)params, vpar, *f, err);
        }
        else {
            double grad_E[wf->getNVP()], dgrad_E[wf->getNVP()];
            if (!f && df) grad(*(struct vmc_workspace *)params, vpar, grad_E, dgrad_E);
            else {
                double en;
                fgrad(*(struct vmc_workspace *)params, vpar, en, err, grad_E, dgrad_E);
                *f = en;
            }
            for (int i=0; i<wf->getNVP(); ++i) gsl_vector_set(df, i, grad_E[i]);
        }
    }

    double vmc_cost(const gsl_vector *v, void *params)
    {
        double result;
        vmc_cost_grad(v, params, &result, NULL);
        return result;
    }

    void vmc_grad (const gsl_vector *v, void *params, gsl_vector *df)
    {
        vmc_cost_grad(v, params, NULL, df);
    }
};

void StochasticReconfigurationOptimization::optimizeWF()
{
    using namespace sropt_details;

    const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_steepest_descent;
    gsl_multimin_fdfminimizer *s = NULL;
    gsl_vector *x;
    gsl_multimin_function_fdf target_func;

    size_t iter = 0;
    int status;

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

    // Initialize method and iterate
    target_func.n = npar;
    target_func.f = vmc_cost;
    target_func.df = vmc_grad;
    target_func.fdf = vmc_cost_grad;
    target_func.params = &w;

    s = gsl_multimin_fdfminimizer_alloc(T, npar);
    gsl_multimin_fdfminimizer_set(s, &target_func, x, 0.5, 0.1);

    do
        {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate(s);

            if (status)
                break;

            status = gsl_multimin_test_gradient(s->gradient, 1e-3);

            if (status == GSL_SUCCESS)
                {
                    printf ("converged to minimum at:\n");
                }

            printf ("%5zu f() = %7.3f, ", iter, s->f);
            for (int i=0; i<npar; ++i) {
                printf("grad %3d = %.3f\n", i, gsl_vector_get(s->gradient, i));
            }
        }
    while (status == GSL_CONTINUE && iter < 1000);

    gsl_vector_free(x);
    gsl_multimin_fdfminimizer_free (s);

    /*// create targetfunction
      StochasticReconfigurationTargetFunction * targetf = new StochasticReconfigurationTargetFunction(_wf, _H, _Nmc, getMCI());
      // declare the Dynamic Descent object
      DynamicDescent * ddesc = new DynamicDescent(targetf);
      // allocate an array that will contain the wave function variational parameters
      double * wfpar = new double[_wf->getNVP()];
      // get the variational parameters
      _wf->getVP(wfpar);
      // set the actual variational parameters as starting point for the Conjugate Gradient algorithm
      ddesc->setX(wfpar);
      // find the optimal parameters by minimizing the energy with the Conjugate Gradient algorithm
      ddesc->findMin();
      // set the found parameters in the wave function
      ddesc->getX(wfpar);
      _wf->setVP(wfpar);
      // free memory
      delete[] wfpar;
      delete ddesc;
      delete targetf;*/
    _Nmc = _Nmc;
}
