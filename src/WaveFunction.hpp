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
    virtual void setVP(const double *vp) = 0;    // --- MUST BE IMPLEMENTED
    virtual void getVP(double *vp) = 0;    // --- MUST BE IMPLEMENTED


    // --- methods herited from MCISamplingFunctionInterface
    //   // wave function values that will be used to compute the acceptance  --- MUST BE IMPLEMENTED
    //   virtual void samplingFunction(const double * in, double * out) = 0;
    //   // MCI acceptance starting from the new and old sampling functions  --- MUST BE IMPLEMENTED
    //   virtual double getAcceptance() = 0;


    // --- computation of the derivatives
    // When called, this method computes all the internal values, such as the derivatives,
    // and stored them internally, ready to be accessed with the getters methods. 
    // It requires the positions as input
    virtual void computeAllInternalValues(const double *) = 0;    // --- MUST BE IMPLEMENTED


    // --- getters for the derivatives
    // first derivative divided by the wf  --- MUST BE IMPLEMENTED
    virtual double getD1LogWF(const int &) = 0; // 0 <= i <= _ndim
    // second derivative divided by the wf  --- MUST BE IMPLEMENTED
    virtual double getD2LogWF(const int &) = 0;  // 0 <= i <= _ndim
    // variational derivative divided by the wf  --- MUST BE IMPLEMENTED
    virtual double getVD1LogWF(const int &) = 0;  // 0 <= i <= _npv
};


#endif
