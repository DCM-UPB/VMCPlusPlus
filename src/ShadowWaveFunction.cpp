#include "ShadowWaveFunction.hpp"

#include <stdexcept>
#include <cmath>


void ShadowWaveFunction::computeAllDerivatives(const double *x){
    /*
    It assumes that, since this function is called on acceptance, the x here is the same as
    the one passed to samplingFunction() last time that it was called.
    Therefore, we can/should use the same shadows coordinates sampled in samplinfFunction().
    */
    const double div2NumSWFSampling = 1./(2.*_num_swf_sampling);
    const double two_div_tau = 2./_tau;
    const double invtau = 1./_tau;

    double * d1 = _getD1DivByWF();
    for (int i=0; i<getTotalNDim(); ++i){
        d1[i] = 0.;
        for (int isampling=0; isampling<_num_swf_sampling; ++isampling){
            d1[i] -= two_div_tau * (x[i] - _s1[isampling][i]);
            d1[i] -= two_div_tau * (x[i] - _s2[isampling][i]);
        }
        d1[i] *= div2NumSWFSampling;
    }

    double * d2 = _getD2DivByWF();
    for (int i=0; i<getTotalNDim(); ++i){
        d2[i] = 0.;
        for (int isampling=0; isampling<_num_swf_sampling; ++isampling){
            d2[i] += pow(x[i] - _s1[isampling][i], 2);
            d2[i] += pow(x[i] - _s2[isampling][i], 2);
        }
        d2[i] *= pow(two_div_tau, 2);
        d2[i] *= div2NumSWFSampling;
        d2[i] -= two_div_tau;
    }

    if (hasVD1()){
        double * vd1 = _getVD1DivByWF();
        // kernel component
        vd1[0] = 0.;
        for (int i=0; i<getTotalNDim(); ++i){
            for (int isampling=0; isampling<_num_swf_sampling; ++isampling){
                vd1[0] += pow(x[i] - _s1[isampling][i], 2);
                vd1[0] += pow(x[i] - _s2[isampling][i], 2);
            }
        }
        vd1[0] = vd1[0]*div2NumSWFSampling * pow(invtau, 2);
        // pure shadow components
        for (int ivp=1; ivp<getNVP(); ++ivp){
            vd1[ivp] = 0.;   // initialize to zero
        }
        int contvp = 1;
        for (PureShadowWaveFunction * pswf : _pswfs){
            for (int isampling=0; isampling<_num_swf_sampling; ++isampling){
                pswf->computeAllDerivatives(_s1[isampling]);
                for (int ivp=0; ivp<pswf->getNVP(); ++ivp){
                    vd1[ivp+contvp] += pswf->getVD1DivByWF(ivp);
                }
                pswf->computeAllDerivatives(_s2[isampling]);
                for (int ivp=0; ivp<pswf->getNVP(); ++ivp){
                    vd1[ivp+contvp] += pswf->getVD1DivByWF(ivp);
                }
            }
            contvp += pswf->getNVP();
        }
        for (int ivp=1; ivp<getNVP(); ++ivp){
            vd1[ivp] *= div2NumSWFSampling;   // divide by the number of samplings
        }
    }

    if (hasD1VD1()){
        double ** d1vd1 = _getD1VD1DivByWF();
        for (int i=0; i<getTotalNDim(); ++i){
            for (int ivp=0; ivp<getNVP(); ++ivp){
                d1vd1[i][ivp] = getD1DivByWF(i) * getVD1DivByWF(ivp);
            }
            d1vd1[i][0] -= getD1DivByWF(0)*invtau;
        }
    }

    if (hasD2VD1()){
        double ** d2vd1 = _getD2VD1DivByWF();
        for (int i=0; i<getTotalNDim(); ++i){
            for (int ivp=0; ivp<getNVP(); ++ivp){
                d2vd1[i][ivp] = getD2DivByWF(i) * getVD1DivByWF(ivp);
            }
            // only Tau term
            d2vd1[i][0] += - 2.*pow(getD1DivByWF(i), 2)*invtau + two_div_tau*invtau;
            // SLOWER AND SAME RESULT
            // double xisi2 = 0.;
            // double xs = 0.;
            // double xisi2xs2 = 0.;
            // for (int isampling=0; isampling<_num_swf_sampling; ++isampling){
            //     const double xisi12 = pow(x[i] - _s1[isampling][i], 2);
            //     const double xisi22 = pow(x[i] - _s2[isampling][i], 2);
            //     xisi2 += xisi12;
            //     xisi2 += xisi22;
            //     double xs1 = 0.;
            //     double xs2 = 0.;
            //     for (int idim=0; idim<getTotalNDim(); ++idim){
            //         xs1 += pow(x[idim] - _s1[isampling][idim], 2);
            //         xs2 += pow(x[idim] - _s2[isampling][idim], 2);
            //     }
            //     xs += xs1;
            //     xs += xs2;
            //     xisi2xs2 += xisi12 * xs1;
            //     xisi2xs2 += xisi22 * xs2;
            // }
            // xisi2 *= div2NumSWFSampling;
            // xs *= div2NumSWFSampling;
            // xisi2xs2 *= div2NumSWFSampling;
            //
            // d2vd1[i][0] += pow(two_div_tau, 2)*pow(invtau, 2)*xisi2xs2 - two_div_tau*xs*pow(invtau, 2) - 2.*pow(two_div_tau, 2)*xisi2 + two_div_tau*invtau;
        }
    }
}



void ShadowWaveFunction::samplingFunction(const double * x, double * proto){
    using namespace std;

    // sum of the pure shadow wave functions
    double sum_wf_s1 = 0.;
    double sum_wf_s2 = 0.;

    // normal distribution
    normal_distribution<double> norm;
    const double sigma = sqrt(0.5*_tau);

    // sample and store the shadows (will be used for computing the derivatives)
    for (int isampling=0; isampling<_num_swf_sampling; ++isampling){
        for (int i=0; i<getTotalNDim(); ++i){
            norm = normal_distribution<double>(x[i], sigma);
            _s1[isampling][i] = norm(_rgen);
            _s2[isampling][i] = norm(_rgen);
        }
    }

    if (_pswfs.size() < 1){
        // if there are no pure shadow wf, set the sums to 1
        sum_wf_s1 = 1.;
        sum_wf_s2 = 1.;
    } else {
        // compute the pure shadow wf values and sum them up
        for (int isampling=0; isampling<_num_swf_sampling; ++isampling){
            for (PureShadowWaveFunction * pswf : _pswfs){
                sum_wf_s1 += pswf->value(_s1[isampling]);
                sum_wf_s2 += pswf->value(_s2[isampling]);
            }
        }
    }

    proto[0] = sum_wf_s1 * sum_wf_s2;
}


double ShadowWaveFunction::getAcceptance(const double * protoold, const double * protonew){
    return protonew[0]/protoold[0];
}



void ShadowWaveFunction::getVP(double *vp){
    vp[0] = _tau;
    int contvp = 1;
    for (PureShadowWaveFunction * pswf : _pswfs){
        pswf->getVP(vp+contvp);
        contvp += pswf->getNVP();
    }
}


void ShadowWaveFunction::setVP(const double *vp){
    _tau = vp[0];
    int contvp = 1;
    for (PureShadowWaveFunction * pswf : _pswfs){
        pswf->setVP(vp+contvp);
        contvp += pswf->getNVP();
    }
}



void ShadowWaveFunction::addPureShadowWaveFunction(PureShadowWaveFunction * pswf){
    if (pswf->getNSpaceDim() != getNSpaceDim()){
        throw std::invalid_argument( "Provided pure shadow wf's getNSpaceDim() is not valid" );
    }
    if (pswf->getNPart() != getNPart()){
        throw std::invalid_argument( "Provided pure shadow wf's getNPart() is not valid" );
    }
    if (pswf->hasVD1() != hasVD1()){
        throw std::invalid_argument( "Provided pure shadow wf's hasVD1() is not valid" );
    }
    _pswfs.push_back(pswf);
    setNVP( getNVP() + pswf->getNVP() );
}





ShadowWaveFunction::ShadowWaveFunction(const double &tau, const int &num_swf_sampling, const int &nspacedim, const int &npart, bool flag_vd1, bool flag_d1vd1, bool flag_d2vd1):
WaveFunction(nspacedim, npart, 1, 1, flag_vd1, flag_d1vd1, flag_d2vd1){
    if (tau <= 0.){
        throw std::invalid_argument( "The Shadow Wave Function parameter tau must be strictly greater than zero" );
    }
    if (num_swf_sampling <= 0){
        throw std::invalid_argument( "The Shadow Wave Function parameter num_swf_sampling must be strictly greater than zero" );
    }
    if (hasD1VD1() && !hasVD1()){
        throw std::invalid_argument( "ShadowWaveFunction derivative d1vd1 requires vd1" );
    }
    if (hasD2VD1() && !hasVD1()){
        throw std::invalid_argument( "ShadowWaveFunction derivative d2vd1 requires vd1 and d1vd1" );
    }

    _tau = tau;
    _num_swf_sampling = num_swf_sampling;
    _s1 = new double*[_num_swf_sampling];
    for (int i=0; i<_num_swf_sampling; ++i) _s1[i] = new double[getTotalNDim()];
    _s2 = new double*[_num_swf_sampling];
    for (int i=0; i<_num_swf_sampling; ++i) _s2[i] = new double[getTotalNDim()];
    _rgen = std::mt19937_64(_rdev());
}



ShadowWaveFunction::~ShadowWaveFunction(){
    for (int i=0; i<_num_swf_sampling; ++i) delete[] _s1[i];
    delete[] _s1;
    for (int i=0; i<_num_swf_sampling; ++i) delete[] _s2[i];
    delete[] _s2;
    _pswfs.clear();
}
