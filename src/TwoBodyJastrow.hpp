#ifndef TWO_BODY_JASTROW
#define TWO_BODY_JASTROW

#include "WaveFunction.hpp"
#include "TwoBodyPseudoPotential.hpp"
#include "ParticleArrayHelper.hpp"

#include <stdexcept>


/*
TwoBodyJastrow is a virtual class (or interface) for any 2-body Jastrow of the form:

    J(R) = exp( sum_ij u(r_ij) )

where R contains all the particle coordinates, sum_ij is a sum over all the particle pairs (i, j) and r_ij is the distance between the particles i and j.
u is a 2-body pseudopotential, and must be implemented as TwoBodyPseudoPotential.
*/


class TwoBodyJastrow: public WaveFunction{

private:
    TwoBodyPseudoPotential * _u2;
    ParticleArrayHelper * _pah;

public:
    TwoBodyJastrow(const int &npart, TwoBodyPseudoPotential * u2):
    WaveFunction(u2->getNSpaceDim(), npart, 1, u2->getNVP(), u2->hasVD1(), u2->hasD1VD1(), u2->hasD2VD1()){
        _u2 = u2;
        _pah = new ParticleArrayHelper(u2->getNSpaceDim());
        if (hasD1VD1() && !hasVD1()){
            throw std::invalid_argument( "TwoBodyJastrow derivative d1vd1 requires vd1" );
        }
        if (hasD2VD1() && !(hasVD1() && hasD1VD1())){
            throw std::invalid_argument( "TwoBodyJastrow derivative d2vd1 requires vd1 and d1vd1" );
        }
    }
    virtual ~TwoBodyJastrow(){
        delete _pah;
    }



    void setVP(const double *vp){_u2->setVP(vp);}
    void getVP(double *vp){_u2->getVP(vp);}



    void samplingFunction(const double * x, double * protov);

    double getAcceptance(const double * protoold, const double * protonew);

    void computeAllDerivatives(const double *x);

    double computeWFValue(const double * protovalues);
};






























#endif
