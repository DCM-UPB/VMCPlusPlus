#ifndef WAVE_FUNCTION
#define WAVE_FUNCTION

#include "MCISamplingFunctionInterface.hpp"


class WaveFunction: public MCISamplingFunctionInterface
{
protected:
    int _nvp;  //number of variational parameters involved
    int _npart;
    int _nspacedim;

public:
    WaveFunction(const int &nspacedim, const int &npart, const int &ncomp, const int &nvp): MCISamplingFunctionInterface(nspacedim*npart,ncomp)
    {
        _nspacedim=nspacedim;
        _npart=npart;
        _nvp=nvp;
    }
    virtual ~WaveFunction(){}

    int getNSpaceDim(){return _nspacedim;}
    int getNPart(){return _npart;}
    int getNVP(){return _nvp;}

    // --- interface for manipulating the variational parameters
    virtual void setVP(const double *vp) = 0;
    virtual void getVP(double *vp) = 0;

    // --- methods herited from MCISamplingFunctionInterface
    //   // wave function values that will be used to compute the acceptance  --- MUST BE IMPLEMENTED
    //   virtual void samplingFunction(const double * in, double * out) = 0;
    //   // MCI acceptance starting from the new and old sampling functions  --- MUST BE IMPLEMENTED
    //   virtual double getAcceptance() = 0;

    //// value of the wf itself  --- MUST BE IMPLEMENTED
    //virtual double wf(const double *) = 0;
    // first derivative divided by the wf  --- MUST BE IMPLEMENTED
    virtual double d1(const int &, const double * ) = 0; // 0 <= i <= _ndim
    // second derivative divided by the wf  --- MUST BE IMPLEMENTED
    virtual double d2(const int &, const double *) = 0;  // 0 <= i <= _ndim
    // variational derivative divided by the wf  --- MUST BE IMPLEMENTED
    virtual double vd1(const int &, const double *) = 0;  // 0 <= i <= _npv
};


#endif
