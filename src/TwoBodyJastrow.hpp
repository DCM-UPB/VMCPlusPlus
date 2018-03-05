#ifndef TWO_BODY_JASTROW
#define TWO_BODY_JASTROW

#include "TwoBodyPseudoPotential.hpp"
#include "ParticlesPosition.hpp"


/*
TwoBodyJastrow is a virtual class (or interface) for any 2-body Jastrow of the form:

    J(R) = exp( - sum_ij u(r_ij) )

where R contains all the particle coordinates, sum_ij is a sum over all the particle pairs (i, j) and r_ij is the distance between the particles i and j.
u is a 2-body pseudopotential, and must be implemented as TwoBodyPseudoPotential.
*/


class TwoBodyJastrow: public WaveFunction{

private:
    TwoBodyPseudoPotential * _u2;
    ParticlesPosition * _ppos;

public:
    TwoBodyJastrow(TwoBodyPseudoPotential * u2, bool flag_vd1=true, bool flag_d1vd1=true, bool flag_d2vd1=true):
    WaveFunction(u2->getNSpaceDim(), u2->getNPart(), 1, u2->getNVP(), flag_vd1, flag_d1vd1, flag_d2vd1){
        _u2 = u2;
        _ppos = new ParticlesPosition(u2->getNSpaceDim());
    }
    virtual ~TwoBodyJastrow();



    void setVP(const double *vp){u2->setVP(vp);}
    void getVP(double *vp){u2->getVP(vp);}



    void samplingFunction(const double * x, double * protov);

    double getAcceptance();

    void computeAllDerivatives(const double *x);

};






























#endif
