#include "vmc/EnergyGradientTargetFunction.hpp"

#include "vmc/EnergyGradientMCObservable.hpp"
#include "vmc/MPIVMC.hpp"

namespace vmc
{

nfm::NoisyValue EnergyGradientTargetFunction::f(const std::vector<double> &vp)
{
    // set the variational parameters given as input
    _vmc.setVP(vp.data());
    // perform the integral and store the values
    double obs[4];
    double dobs[4];
    _vmc.computeEnergy(_E_Nmc, obs, dobs, true, true);
    MPIVMC::Integrate(_vmc.getMCI(), _E_Nmc, obs, dobs, true, true);
    nfm::NoisyValue f{obs[0], dobs[0]};

    if (_lambda_reg > 0.) { // compute the regularization term
        const double norm = std::inner_product(vp.begin(), vp.end(), vp.begin(), 0.);
        f.val += _lambda_reg*norm/_vmc.getNVP();
    }

    return f;
}


void EnergyGradientTargetFunction::grad(const std::vector<double> &vp, nfm::NoisyGradient &grad)
{
    fgrad(vp, grad);
}


nfm::NoisyValue EnergyGradientTargetFunction::fgrad(const std::vector<double> &vp, nfm::NoisyGradient &grad)
{
    const int nvp = _vmc.getNVP();

    // set the variational parameters given as input
    _vmc.setVP(vp.data());
    // add gradient obs to MCI
    const int blocksize = this->hasGradErr() ? 1 : 0;
    _vmc.getMCI().addObservable(EnergyGradientMCObservable(_vmc.getNTotalDim(), nvp), blocksize, 1, false, blocksize > 0); // skipping equlibiration for gradients
    // perform the integral and store the values
    double obs[4 + 2*nvp];
    double dobs[4 + 2*nvp];
    _vmc.computeEnergy(_grad_E_Nmc, obs, dobs, true, true);
    _vmc.getMCI().popObservable(); // remove the gradient obs (it will be deleted)
    // create pointers for ease of use and readability
    const double * const H = obs;
    const double * const dH = dobs;
    const double * const Oi = obs + 4;
    const double * const dOi = dobs + 4;
    const double * const HOi = obs + 4 + nvp;
    const double * const dHOi = dobs + 4 + nvp;
    nfm::NoisyValue f{H[0], dH[0]};
    // compute direction (or gradient) to follow
    for (int i = 0; i < nvp; ++i) {
        grad.val[i] = -2.*(HOi[i] - H[0]*Oi[i]);
    }
    if (this->hasGradErr()) {
        for (int i = 0; i < nvp; ++i) { // error propagation
            grad.err[i] = 2.*sqrt(dHOi[i]*dHOi[i] + dH[0]*dH[0]*Oi[i]*Oi[i] + H[0]*H[0]*dOi[i]*dOi[i]);
        }
    }

    if (_lambda_reg > 0.) { // compute the regularization terms
        const double norm = std::inner_product(vp.begin(), vp.end(), vp.begin(), 0.);
        const double fac = _lambda_reg/nvp;
        f.val += fac*norm;
        for (int i = 0; i < nvp; ++i) { grad.val[i] += 2.*fac*vp[i]; }
    }

    return f;
}
} // namespace vmc