#include "FunctionWithVariationalParameters.hpp"

#include <iostream>

using namespace std;



int FunctionWithVariationalParameters::getNVP(){
    return _nvp;
}


void FunctionWithVariationalParameters::setNVP(const int &nvp){
    _nvp = nvp;
    delete[] _vp;
    _vp = new double[_nvp];
}


void FunctionWithVariationalParameters::getVP(double *vp){
    for (int i=0; i<getNVP(); ++i){
        vp[i] = _vp[i];
    }
}


double FunctionWithVariationalParameters::getVP(const int &i){
    return _vp[i];
}


void FunctionWithVariationalParameters::setVP(const double *vp){
    for (int i=0; i<getNVP(); ++i){
        _vp[i] = vp[i];
        actAfterVPChange(i, vp[i]);
    }
}


void FunctionWithVariationalParameters::setVP(const int &i, const double &vp){
    _vp[i] = vp;
    actAfterVPChange(i, vp);
}



FunctionWithVariationalParameters::FunctionWithVariationalParameters(const int &nvp){
    _nvp = nvp;
    _vp = new double[nvp];
}


FunctionWithVariationalParameters::~FunctionWithVariationalParameters(){
    _nvp = 0;
    delete[] _vp;
    _vp = 0;
}
