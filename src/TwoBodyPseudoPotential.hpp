#ifndef TWO_BODY_PSEUDO_POTENTIAL
#define  TWO_BODY_PSEUDO_POTENTIAL


#include "ParticleDistance.hpp"



class TwoBodyPseudoPotential{

private:
    ParticleDistance * _dist;

    double * _foo;
    double * _foo2;
    double * _vfoo;
    double * _vfoo2;

public:
    TwoBodyPseudoPotential(ParticleDistance * dist){
        _dist = dist;
        _foo = new double[2*_dist->getNSpaceDim()];
        _foo2 = 0;
        _vfoo = 0;
        _vfoo2 = 0;
    }
    virtual ~TwoBodyPseudoPotential(){
        delete[] _foo;
        if (_foo2) delete[] _foo2;
        if (_vfoo) delete[] _vfoo;
        if (_vfoo2) delete[] _vfoo2;
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
    virtual void urD2VD1(const double &dist, double * d1vd1) = 0;



    // -- derivatives
    double u(const double * r1, const double * r2);
    void d1(const double * r1, const double * r2, double * deriv1);
    void d2(const double * r1, const double * r2, double * deriv2);
    void vd1(const double * r1, const double * r2, double * varderiv1);
    void d1vd1(const double * r1, const double * r2, double ** deriv1varderiv1);
    void d2vd1(const double * r1, const double * r2, double ** deriv2varderiv1);

};


#endif
