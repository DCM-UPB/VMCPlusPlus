#ifndef TWO_BODY_PSEUDO_POTENTIAL
#define  TWO_BODY_PSEUDO_POTENTIAL


#include "ParticleDistance.hpp"



class TwoBodyPseudoPotential{

private:
    ParticleDistance * _dist;
    int _npart;

    double * _foo;

public:
    TwoBodyPseudoPotential(const int &npart, ParticleDistance * dist){
        _npart = npart;
        _dist = dist;
        _foo = new double[_dist->getNSpaceDim()];
    }
    virtual ~TwoBodyPseudoPotential(){
        delete[] _foo;
    }

    int getNSpaceDim(){return _dist->getNSpaceDim();}
    int getNPart(){return _npart;}


    virtual int getNVP() = 0;
    virtual void setVP(const double *vp) = 0;
    virtual void getVP(double *vp) = 0;

    virtual double ur(const double &dist) = 0;
    virtual double urD1(const double &dist) = 0;
    virtual double urD2(const double &dist) = 0;


    double u(const double * r1, const double * r2);
    void d1(const double * r1, const double * r2, double * deriv1);
    void d2(const double * r1, const double * r2, double * deriv2);

    virtual double vd1() = 0;
    virtual double d1vd1() = 0;
    virtual double d2vd1() = 0;

};


#endif
