#include "WaveFunction.hpp"




void WaveFunction::setNVP(const int &nvp){
    _nvp=nvp;
    _allocateVariationalDerivativesMemory();
}



void WaveFunction::callBackFunction(const double *x, const bool flag_observation){
    if (flag_observation){
        computeAllDerivatives(x);
    }
}


void WaveFunction::_allocateVariationalDerivativesMemory(){
    if (hasVD1()){
        if (_vd1_divbywf != 0){
            delete[] _vd1_divbywf;
        }
        _vd1_divbywf = new double[getNVP()];
        for (int i=0; i<getNVP(); ++i) _vd1_divbywf[i] = 0.;
    }
    if (hasD1VD1()){
        if (_d1vd1_divbywf != 0){
            for (int i=0; i<getTotalNDim(); ++i){
                delete[] _d1vd1_divbywf[i];
            }
            delete[] _d1vd1_divbywf;
        }
        _d1vd1_divbywf = new double*[getTotalNDim()];
        for (int i=0; i<getTotalNDim(); ++i){
            _d1vd1_divbywf[i] = new double[getNVP()];
            for (int j=0; j<getNVP(); ++j) _d1vd1_divbywf[i][j] = 0.;
        }
    }
    if (hasD1VD1()){
        if (_d2vd1_divbywf != 0){
            for (int i=0; i<getTotalNDim(); ++i){
                delete[] _d2vd1_divbywf[i];
            }
            delete[] _d2vd1_divbywf;
        }
        _d2vd1_divbywf = new double*[getTotalNDim()];
        for (int i=0; i<getTotalNDim(); ++i){
            _d2vd1_divbywf[i] = new double[getNVP()];
            for (int j=0; j<getNVP(); ++j) _d2vd1_divbywf[i][j] = 0.;
        }
    }
}


WaveFunction::WaveFunction(const int &nspacedim, const int &npart, const int &ncomp, const int &nvp, bool flag_vd1, bool flag_d1vd1, bool flag_d2vd1):
MCISamplingFunctionInterface(nspacedim*npart, ncomp),
MCICallBackOnAcceptanceInterface(nspacedim*npart){
    _nspacedim=nspacedim;
    _npart=npart;
    _nvp=nvp;

    _flag_vd1 = flag_vd1;
    _flag_d1vd1 = flag_d1vd1;
    _flag_d2vd1 = flag_d2vd1;

    _d1_divbywf = new double[nspacedim*npart];
    _d2_divbywf = new double[nspacedim*npart];

    for (int i=0; i<nspacedim*npart; ++i){
        _d1_divbywf[i] = 0.;
        _d2_divbywf[i] = 0.;
    }

    _vd1_divbywf = 0;
    _d1vd1_divbywf = 0;
    _d2vd1_divbywf = 0;
    _allocateVariationalDerivativesMemory();
}


WaveFunction::~WaveFunction(){
    delete[] _d1_divbywf; _d1_divbywf = 0;
    delete[] _d2_divbywf; _d2_divbywf = 0;
    if (_vd1_divbywf != 0){
        delete[] _vd1_divbywf; _vd1_divbywf = 0;
    }
    if (_d1vd1_divbywf != 0){
        for (int i=0; i<getTotalNDim(); ++i){
            delete[] _d1vd1_divbywf[i];
        }
        delete[] _d1vd1_divbywf; _d1vd1_divbywf = 0;
    }
    if (_d2vd1_divbywf != 0){
        for (int i=0; i<getTotalNDim(); ++i){
            delete[] _d2vd1_divbywf[i];
        }
        delete[] _d2vd1_divbywf; _d2vd1_divbywf = 0;
    }
}
