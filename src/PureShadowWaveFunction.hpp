#ifndef PURE_SHADOW_WAVE_FUNCTION
#define PURE_SHADOW_WAVE_FUNCTION

#include "FunctionWithNSpaceDimAndNPart.hpp"
#include "FunctionWithVariationalParameters.hpp"
#include "FunctionWithDerivatives.hpp"

/*
IMPLEMENTATIONS OF THIS INTERFACE MUST INCLUDE:

    - void getVP(double *vp)
            get the variational parameters

    - void setVP(const double *vp)
            set the variational parameters

    - double value(const double *s)
            return the value of the Pure Shadow Wave Function as function of the shadow position s

    - void computeAllDerivatives(const double *x)
            Only for computing (if required) the variational derivatives VD1. Use the setters for derivatives values (setVD1DivByWF)


*/



class PureShadowWaveFunction: public FunctionWithNSpaceDimAndNPart, public FunctionWithVariationalParameters, public FunctionWithDerivatives{
public:
    PureShadowWaveFunction(const int &nspacedim, const int &npart, const int &nvp, bool flag_vd1):
    FunctionWithNSpaceDimAndNPart(nspacedim, npart),
    FunctionWithVariationalParameters(nvp),
    FunctionWithDerivatives(nspacedim*npart, nvp, false, false, flag_vd1, false, false){

    }

    ~PureShadowWaveFunction(){

    }


    // --- value of the the Pure Shadow Wave Function Component
    virtual double value(const double *s) = 0;

};


#endif
