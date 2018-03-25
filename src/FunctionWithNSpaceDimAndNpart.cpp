#include "FunctionWithNSpaceDimAndNPart.hpp"



int FunctionWithNSpaceDimAndNPart::getNSpaceDim(){
    return _nspacedim;
}



int FunctionWithNSpaceDimAndNPart::getTotalNDim(){
    return _nspacedim * _npart;
}



int FunctionWithNSpaceDimAndNPart::getNPart(){
    return _npart;
}



FunctionWithNSpaceDimAndNPart::FunctionWithNSpaceDimAndNPart(const int &nspacedim, const int &npart){
    _nspacedim = nspacedim;
    _npart = npart;
}



FunctionWithNSpaceDimAndNPart::~FunctionWithNSpaceDimAndNPart(){
    _nspacedim = 0;
    _npart = 0;
}
