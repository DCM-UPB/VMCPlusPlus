#ifndef TWO_BODY_JASTROW
#define TWO_BODY_JASTROW

#include "WaveFunction.hpp"
#include "TwoBodyPseudoPotential.hpp"
#include "ParticleArrayHelper.hpp"

#include <stdexcept>
#include <iostream>
using namespace std;


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
    TwoBodyJastrow(const int &npart, TwoBodyPseudoPotential * u2);

    virtual ~TwoBodyJastrow(){
        delete _pah;
    }

    void getVP(double * vp){_u2->getVP(vp);}
    void setVP(const double * vp){_u2->setVP(vp);}

    void samplingFunction(const double * x, double * protov);

    double getAcceptance(const double * protoold, const double * protonew);

    void computeAllDerivatives(const double *x);

};






























#endif
