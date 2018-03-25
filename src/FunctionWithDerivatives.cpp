#include "FunctionWithDerivatives.hpp"




void FunctionWithDerivatives::_allocateDerivativesMemory(const int &totalndim, const int &nvp){

    if (hasD1() && (totalndim != _totalndim)){
        if ( _d1_divbywf != 0){
            delete[] _d1_divbywf;
        }
        _d1_divbywf = new double[totalndim];
    }

    if (hasD2() && (totalndim != _totalndim)){
        if ( _d2_divbywf != 0){
            delete[] _d2_divbywf;
        }
        _d2_divbywf = new double[totalndim];
    }

    if (hasVD1() && (nvp != _nvp)){
        if ( _vd1_divbywf != 0 ){
            delete[] _vd1_divbywf;
        }
        _vd1_divbywf = new double[nvp];
    }

    if (hasD1VD1() && ((totalndim != _totalndim) || (nvp != _nvp))){
        if ( _d1vd1_divbywf != 0 ){
            for (int i=0; i<_totalndim; ++i){
                delete[] _d1vd1_divbywf[i];
            }
            delete[] _d1vd1_divbywf;
        }
        _d1vd1_divbywf = new double*[totalndim];
        for (int i=0; i<totalndim; ++i){
            _d1vd1_divbywf[i] = new double[nvp];
        }
    }

    if (hasD2VD1() && ((totalndim != _totalndim) || (nvp != _nvp))){
        if ( _d2vd1_divbywf != 0 ){
            for (int i=0; i<_totalndim; ++i){
                delete[] _d2vd1_divbywf[i];
            }
            delete[] _d2vd1_divbywf;
        }
        _d2vd1_divbywf = new double*[totalndim];
        for (int i=0; i<totalndim; ++i){
            _d2vd1_divbywf[i] = new double[nvp];
        }
    }

    _totalndim = totalndim;
    _nvp = nvp;
}



FunctionWithDerivatives::FunctionWithDerivatives(const int &totalndim, const int &nvp, bool flag_d1, bool flag_d2, bool flag_vd1, bool flag_d1vd1, bool flag_d2vd1){
    _flag_d1 = flag_d1;
    _flag_d2 = flag_d2;
    _flag_vd1 = flag_vd1;
    _flag_d1vd1 = flag_d1vd1;
    _flag_d2vd1 = flag_d2vd1;

    _d1_divbywf = 0;
    _d2_divbywf = 0;
    _vd1_divbywf = 0;
    _d1vd1_divbywf = 0;
    _d2vd1_divbywf = 0;
    _totalndim = 0;
    _nvp = 0;

    _allocateDerivativesMemory(totalndim, nvp);
}



FunctionWithDerivatives::~FunctionWithDerivatives(){
    if (hasD1()){
        delete[] _d1_divbywf;
        _d1_divbywf = 0;
    }

    if (hasD2()){
        delete[] _d2_divbywf;
        _d2_divbywf = 0;
    }

    if (hasVD1()){
        delete[] _vd1_divbywf;
        _vd1_divbywf = 0;
    }

    if (hasD1VD1()){
        for (int i=0; i<_totalndim; ++i){
            delete[] _d1vd1_divbywf[i];
        }
        delete[] _d1vd1_divbywf;
        _d1vd1_divbywf = 0;
    }

    if (hasD2VD1()){
        for (int i=0; i<_totalndim; ++i){
            delete[] _d2vd1_divbywf[i];
        }
        delete[] _d2vd1_divbywf;
        _d2vd1_divbywf = 0;
    }
}
