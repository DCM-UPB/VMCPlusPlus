#include "vmc/StochasticReconfigurationTargetFunction.hpp"
#include "vmc/StochasticReconfigurationMCObservable.hpp"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

namespace vmc
{

void add_norm_f(const double * const vp, const int nvp, double &f, const double lambda)
{
    // compute the normalization term
    const double norm = std::inner_product(vp, vp + nvp, vp, 0.);
    f += lambda*norm/nvp;
}

void add_norm_grad(const double * const vp, const int nvp, double * const grad_E, const double lambda)
{
    const double fac = lambda/nvp;
    // add the normalization gradient
    for (int i = 0; i < nvp; ++i) {
        grad_E[i] += 2.*vp[i]*fac;
    }
}

void add_norm_fgrad(const double * const vp, const int nvp, double &f, double * const grad_E, const double lambda)
{
    const double fac = lambda/nvp;

    // compute the normalization term
    const double norm = std::inner_product(vp, vp + nvp, vp, 0.);

    // add the normalization value and gradient
    f += norm*fac;
    for (int i = 0; i < nvp; ++i) {
        grad_E[i] += 2.*vp[i]*fac;
    }
}


void StochasticReconfigurationTargetFunction::_integrate(const double * const vp, double * const obs, double * const dobs, const bool flag_grad, const bool flag_dgrad)
{
    // set the variational parameters given as input
    _vmc.setVP(vp);

    // set up the MC integrator
    if (flag_grad) { // add gradient obs if necessary
        // skip MC error for grad if flag_dgrad is false
        _vmc.getMCI().addObservable(StochasticReconfigurationMCObservable(_vmc.getNTotalDim(), _vmc.getNVP()),
                                    flag_dgrad ? 1 : 0, 1, false, flag_dgrad);
    }

    // perform the integral and store the values (skip extra burning phase on gradient runs (only findMRT2))
    _vmc.computeEnergy(flag_grad ? _grad_E_Nmc : _E_Nmc, obs, dobs, true, !flag_grad);

    // remove gradient obs again
    if (flag_grad) { _vmc.getMCI().popObservable(); }
}

void StochasticReconfigurationTargetFunction::_calcObs(const double * const vp, double &f, double &df, double * const grad_E, double * const dgrad_E)
{
    const auto nvp = static_cast<size_t>(_vmc.getNVP());
    const bool flag_grad = (grad_E != nullptr);
    const bool flag_dgrad = (dgrad_E != nullptr);

    double obs[flag_grad ? 4 + 2*nvp + nvp*nvp : 4];
    double dobs[flag_grad ? 4 + 2*nvp + nvp*nvp : 4];

    _integrate(vp, obs, dobs, flag_grad, flag_dgrad);

    f = obs[0];
    df = dobs[0];

    if (flag_grad) {
        // create pointers for ease of use and readability
        double * const H = obs;
        double * const dH = dobs;
        double * const Oi = obs + 4;
        double * const dOi = dobs + 4;
        double * const HOi = obs + 4 + nvp;
        double * const dHOi = dobs + 4 + nvp;
        double * const OiOj = obs + 4 + 2*nvp;
        double * const dOiOj = dobs + 4 + 2*nvp;


        // --- compute direction (or gradient) to follow
        gsl_matrix * sij = gsl_matrix_alloc(nvp, nvp);
        gsl_matrix * rdsij = flag_dgrad ? gsl_matrix_alloc(nvp, nvp) : nullptr;   // relative error, i.e. error/value
        for (size_t i = 0; i < nvp; ++i) {
            for (size_t j = 0; j < nvp; ++j) {
                gsl_matrix_set(sij, i, j, OiOj[i*nvp + j] - Oi[i]*Oi[j]);
                if (flag_dgrad) {
                    gsl_matrix_set(rdsij, i, j,
                                   (dOiOj[i*nvp + j] + fabs(Oi[i]*Oi[j])*((dOi[i]/Oi[i]) + (dOi[j]/Oi[j])))
                                   /gsl_matrix_get(sij, i, j));
                }
            }
        }
        gsl_vector * fi = gsl_vector_alloc(nvp);
        gsl_vector * rdfi = flag_dgrad ? gsl_vector_alloc(nvp) : nullptr;   // relative error, i.e. error/value
        for (size_t i = 0; i < nvp; ++i) {
            gsl_vector_set(fi, i, H[0]*Oi[i] - HOi[i]);
            if (flag_dgrad) {
                gsl_vector_set(rdfi, i,
                               (fabs(H[0]*Oi[i])*((dH[0]/H[0]) + (dOi[i]/Oi[i])) + dHOi[i])
                               /gsl_vector_get(fi, i));
            }
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
        for (size_t i = 0; i < nvp; ++i) {
            for (size_t j = 0; j < nvp; ++j) {
                gsl_matrix_set(Isij, i, j, 0.);
            }
        }
        for (size_t i = 0; i < nvp; ++i) {
            if (gsl_vector_get(S, i) > SVD_MIN*gsl_vector_get(S, 0)) {
                gsl_matrix_set(Isij, i, i, 1./gsl_vector_get(S, i));
            }
            else {
                gsl_matrix_set(Isij, i, i, 0.);
            }
        }
        gsl_matrix * mm = gsl_matrix_alloc(nvp, nvp);
        gsl_matrix_transpose(V);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Isij, V, 0.0, mm);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, sij, mm, 0.0, Isij);
        // --- finally find the direction to follow
        for (size_t i = 0; i < nvp; ++i) {
            grad_E[i] = 0.;
            if (flag_dgrad) { dgrad_E[i] = 0.; }
            for (size_t k = 0; k < nvp; ++k) {
                double foo = gsl_vector_get(fi, k)*gsl_matrix_get(Isij, k, i);
                grad_E[i] += foo;
                if (flag_dgrad) {
                    dgrad_E[i] += fabs(foo)*(gsl_vector_get(rdfi, k) + gsl_matrix_get(rdsij, k, i));  // not correct, just a rough estimation
                }
            }
        }

        // free resources
        gsl_matrix_free(mm);
        gsl_vector_free(work);
        gsl_vector_free(S);
        gsl_matrix_free(V);
        gsl_matrix_free(Isij);
        if (flag_dgrad) { gsl_vector_free(rdfi); }
        gsl_vector_free(fi);
        if (flag_dgrad) { gsl_matrix_free(rdsij); }
        gsl_matrix_free(sij);
    }
}


nfm::NoisyValue StochasticReconfigurationTargetFunction::f(const std::vector<double> &vp)
{
    nfm::NoisyValue f;
    _calcObs(vp.data(), f.val, f.err);
    if (_lambda_reg > 0) { add_norm_f(vp.data(), _vmc.getNVP(), f.val, _lambda_reg); }
    return f;
}

void StochasticReconfigurationTargetFunction::grad(const std::vector<double> &vp, nfm::NoisyGradient &grad)
{
    double f, df; // dummies
    _calcObs(vp.data(), f, df, grad.val.data(), this->hasGradErr() ? grad.err.data() : nullptr);
    if (_lambda_reg > 0) { add_norm_grad(vp.data(), _vmc.getNVP(), grad.val.data(), _lambda_reg); }
}

nfm::NoisyValue StochasticReconfigurationTargetFunction::fgrad(const std::vector<double> &vp, nfm::NoisyGradient &grad)
{
    nfm::NoisyValue f;
    _calcObs(vp.data(), f.val, f.err, grad.val.data(), this->hasGradErr() ? grad.err.data() : nullptr);
    if (_lambda_reg > 0) { add_norm_fgrad(vp.data(), _vmc.getNVP(), f.val, grad.val.data(), _lambda_reg); }
    return f;
}
} // namespace vmc