#ifndef TWO_BODY_PSEUDO_POTENTIAL
#define  TWO_BODY_PSEUDO_POTENTIAL


#include "ParticleDistance.hpp"



class TwoBodyPseudoPotential{

private:
    ParticleDistance * _dist;

    double * _foo;
    double * _vfoo;

public:
    TwoBodyPseudoPotential(ParticleDistance * dist){
        _dist = dist;
        _foo = new double[2*_dist->getNSpaceDim()];
        _vfoo = 0;
    }
    virtual ~TwoBodyPseudoPotential(){
        delete[] _foo;
        if (_vfoo) delete[] _vfoo;
    }

    int getNSpaceDim(){return _dist->getNSpaceDim();}

    // --- methods that must be implemented
    virtual int getNVP() = 0;
    virtual void setVP(const double *vp) = 0;
    virtual void getVP(double *vp) = 0;

    virtual double ur(const double &dist) = 0;
    virtual double urD1(const double &dist) = 0;
    virtual double urD2(const double &dist) = 0;
    virtual void urVD1(const double &dist, double * vd1) = 0;
    virtual void urD1VD1(const double &dist, double * d1vd1) = 0;



    // -- derivatives
    double u(const double * r1, const double * r2);
    void d1(const double * r1, const double * r2, double * deriv1);
    void d2(const double * r1, const double * r2, double * deriv2);
    void vd1(const double * r1, const double * r2, double * varderiv1);
    void d1vd1(const double * r1, const double * r2, double ** deriv1varderiv1);

    // virtual double d1vd1() = 0;
    // virtual double d2vd1() = 0;

};


#endif
