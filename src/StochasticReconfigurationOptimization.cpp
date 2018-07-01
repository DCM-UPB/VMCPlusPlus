#include "StochasticReconfigurationOptimization.hpp"
#include "StochasticReconfigurationMCObservable.hpp"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

namespace sropt_details {
    void vmc_workspace::initFromOptimizer(StochasticReconfigurationOptimization * wfopt)
    {
        wf = wfopt->getWF();
        H = wfopt->getH();
        mci = wfopt->getMCI();
        Nmc = wfopt->getNmc();
    }

    double vmc_cost(const gsl_vector *v, void *params)
    {
        WaveFunction * const wf = ((struct vmc_workspace *)params)->wf;
        Hamiltonian * const H = ((struct vmc_workspace *)params)->H;
        MCI * const mci = ((struct vmc_workspace *)params)->mci;
        const long Nmc = ((struct vmc_workspace *)params)->Nmc;

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
        mci->integrate(Nmc, energy, d_energy);

        // compute the normalization
        double norm = 0.;
        for (int i=0; i<wf->getNVP(); ++i)
            norm += pow(gsl_vector_get(v, i), 2);
        norm = sqrt(norm)/wf->getNVP();

        // return the cost function
        return energy[0];
    };

    // NoisyFunctionWithGradient implementation
    void f(vmc_workspace &w, const double *vp, double &f, double &df){
        // set the variational parameters given as input
        w.wf->setVP(vp);
        // set up the MC integrator
        w.mci->clearSamplingFunctions(); w.mci->addSamplingFunction(w.wf);
        w.mci->clearObservables(); w.mci->addObservable(w.H);
        // perform the integral and store the values
        double * obs = new double[4];
        double * dobs = new double[4];
        w.mci->integrate(w.Nmc, obs, dobs);
        f = obs[0];
        df = dobs[0];
        // free resources
        delete [] dobs;
        delete [] obs;
    }

    void grad(vmc_workspace &w, const double *vp, double *grad_E, double *dgrad_E){

        // set the variational parameters given as input
        w.wf->setVP(vp);
        // set up the MC integrator
        w.mci->clearSamplingFunctions(); w.mci->addSamplingFunction(w.wf);
        w.mci->clearObservables();
        w.mci->addObservable(w.H);
        StochasticReconfigurationMCObservable * mc_obs = new StochasticReconfigurationMCObservable(w.wf, w.H);
        w.mci->addObservable(mc_obs);
        // perform the integral and store the values
        double * obs = new double[4 + 2*w.wf->getNVP() + w.wf->getNVP()*w.wf->getNVP()];
        double * dobs = new double[4 + 2*w.wf->getNVP() + w.wf->getNVP()*w.wf->getNVP()];
        w.mci->integrate(w.Nmc, obs, dobs);
        // create pointers for ease of use and readability
        double * H = obs;
        double * dH = dobs;
        double * Oi = obs+4;
        double * dOi = dobs+4;
        double * HOi = obs+4+w.wf->getNVP();
        double * dHOi = dobs+4+w.wf->getNVP();
        double * OiOj = obs+4+2*w.wf->getNVP();
        double * dOiOj = dobs+4+2*w.wf->getNVP();

        // --- compute direction (or gradient) to follow
        gsl_matrix * sij = gsl_matrix_alloc(w.wf->getNVP(), w.wf->getNVP());
        gsl_matrix * rdsij = gsl_matrix_alloc(w.wf->getNVP(), w.wf->getNVP());   // relative error, i.e. error/value
        for (int i=0; i<w.wf->getNVP(); ++i){
            for (int j=0; j<w.wf->getNVP(); ++j){
                gsl_matrix_set(sij, i, j, OiOj[i*w.wf->getNVP() + j] - Oi[i] * Oi[j]);
                gsl_matrix_set(rdsij, i, j,
                               (dOiOj[i*w.wf->getNVP() + j] + abs(Oi[i]*Oi[j])*( (dOi[i]/Oi[i]) + (dOi[j]/Oi[j]) ))
                               /  gsl_matrix_get(sij, i, j) );
            }
        }
        gsl_vector * fi = gsl_vector_alloc(w.wf->getNVP());
        gsl_vector * rdfi = gsl_vector_alloc(w.wf->getNVP());   // relative error, i.e. error/value
        for (int i=0; i<w.wf->getNVP(); ++i){
            gsl_vector_set(fi, i, H[0]*Oi[i] - HOi[i]);
            gsl_vector_set(rdfi, i,
                           (abs(H[0]*Oi[i])*( (dH[0]/H[0]) + (dOi[i]/Oi[i]) ) + dHOi[i])
                           / gsl_vector_get(fi, i)  );
        }
        // invert matrix using SVD
        const double SVD_MIN = 1.0e-9;
        // matrix and vectors needed for the SVD
        gsl_matrix * V = gsl_matrix_alloc(w.wf->getNVP(), w.wf->getNVP());
        gsl_vector * S = gsl_vector_alloc(w.wf->getNVP());
        gsl_vector * work = gsl_vector_alloc(w.wf->getNVP());
        // run the Single Value Decomposition
        gsl_linalg_SV_decomp(sij, V, S, work);
        // assemble the inverse matrix
        gsl_matrix * Isij = gsl_matrix_alloc(w.wf->getNVP(), w.wf->getNVP());
        for (int i=0; i<w.wf->getNVP(); ++i){
            for (int j=0; j<w.wf->getNVP(); ++j){
                gsl_matrix_set(Isij, i, j, 0.);
            }
        }
        for (int i=0; i<w.wf->getNVP(); ++i){
            if (gsl_vector_get(S, i) > SVD_MIN * gsl_vector_get(S, 0)){
                gsl_matrix_set(Isij, i, i,   1./gsl_vector_get(S, i)  );
            } else {
                gsl_matrix_set(Isij, i, i,   0.  );
            }
        }
        gsl_matrix * mm = gsl_matrix_alloc(w.wf->getNVP(), w.wf->getNVP());
        gsl_matrix_transpose(V);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Isij, V, 0.0, mm);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, sij, mm, 0.0, Isij);
        // --- finally find the direction to follow
        double foo;
        for (int i=0; i<w.wf->getNVP(); ++i){
            grad_E[i] = 0.;
            dgrad_E[i] = 0.;
            for (int k=0; k<w.wf->getNVP(); ++k){
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

        delete[] dobs;
        delete[] obs;
    }
};

void StochasticReconfigurationOptimization::optimizeWF()
{
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
