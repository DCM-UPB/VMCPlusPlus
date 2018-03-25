#include "WaveFunction.hpp"



void WaveFunction::setNVP(const int &nvp){
    FunctionWithVariationalParameters::setNVP(nvp);
    FunctionWithDerivatives::_allocateDerivativesMemory(getTotalNDim(), getNVP());
}
