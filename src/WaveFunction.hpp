#ifndef WAVE_FUNCTION
#define WAVE_FUNCTION

#include "MCISamplingFunctionInterface.hpp"
#include "MCICallBackOnAcceptanceInterface.hpp"
#include "FunctionWithNSpaceDimAndNPart.hpp"
#include "FunctionWithVariationalParameters.hpp"
#include "FunctionWithDerivatives.hpp"


/*
IMPLEMENTATIONS OF THIS INTERFACE MUST INCLUDE:

    - void samplingFunction(const double * in, double * out)
            heritage from MCISamplingFunctionInterface, uses Psi^2

    - double getAcceptance(const double * protoold, const double * protonew)
            heritage from MCISamplingFunctionInterface

    - void computeAllDerivatives(const double *x)
            use the setters for derivatives values (setD1DivByWF, setD2DivByWF, etc.)

    - void actAfterVPChange(const int &i, const double &vp)


*/


class WaveFunction: public MCISamplingFunctionInterface, public MCICallBackOnAcceptanceInterface, public FunctionWithNSpaceDimAndNPart, public FunctionWithVariationalParameters, public FunctionWithDerivatives{
protected:

public:
    WaveFunction(const int &nspacedim, const int &npart, const int &ncomp, const int &nvp, bool flag_vd1=true, bool flag_d1vd1=true, bool flag_d2vd1=true):
    MCISamplingFunctionInterface(nspacedim*npart, ncomp),
    MCICallBackOnAcceptanceInterface(nspacedim*npart),
    FunctionWithNSpaceDimAndNPart(nspacedim, npart),
    FunctionWithVariationalParameters(nvp),
    FunctionWithDerivatives(nspacedim*npart, nvp, true, true, flag_vd1, flag_d1vd1, flag_d2vd1){}

    virtual ~WaveFunction(){}


    // change the number of variational parameters; it will propagate to FunctionWithVariationalParameters and FunctionWithDerivatives
    void setNVP(const int &nvp){
        FunctionWithVariationalParameters::setNVP(nvp);
        FunctionWithDerivatives::_allocateDerivativesMemory(getTotalNDim(), nvp);
    }


    // --- method inherited from MCICallBackOnAcceptanceInterface, which will simply call computeAllDerivatives
    void callBackFunction(const double *x, const bool flag_observation){
        if (flag_observation){
            computeAllDerivatives(x);
        }
    }



};


#endif
