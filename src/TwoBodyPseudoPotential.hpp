#ifndef TWO_BODY_PSEUDO_POTENTIAL
#define  TWO_BODY_PSEUDO_POTENTIAL


#include "Metric.hpp"



class TwoBodyPseudoPotential{

private:
    Metric * _metric;
    int _ndim2;
    int _nvp;

    bool _flag_vd1;
    bool _flag_d1vd1;
    bool _flag_d2vd1;

    // arrays where the derivatives will be stored
    double * _d1;
    double * _d2;
    double * _vd1;
    double ** _d1vd1;
    double ** _d2vd1;

    // arrays used for computations
    double * _foo;
    double * _foo2;
    double * _vfoo;
    double * _vfoo1;
    double * _vfoo2;

public:
    TwoBodyPseudoPotential(Metric * metric, const int &nvp, bool flag_vd1=true, bool flag_d1vd1=true, bool flag_d2vd1=true);
    virtual ~TwoBodyPseudoPotential();

    int getNSpaceDim(){return _metric->getNSpaceDim();}
    int getNVP(){return _nvp;}



    double u(const double * r1, const double * r2);


    void computeAllDerivatives(const double * r1, const double * r2);




    // -- derivatives
    double getD1(const int &id1){return _d1[id1];}
    double getD2(const int &id2){return _d2[id2];}
    double getVD1(const int &ivd1){return _vd1[ivd1];}
    double getD1VD1(const int &id1, const int &ivd1){return _d1vd1[id1][ivd1];}
    double getD2VD1(const int &id2, const int &ivd1){return _d2vd1[id2][ivd1];}




    // --- methods that must be implemented
    virtual void setVP(const double *vp) = 0;
    virtual void getVP(double *vp) = 0;

    virtual double ur(const double &dist) = 0;
    virtual double urD1(const double &dist) = 0;
    virtual double urD2(const double &dist) = 0;
    virtual void urVD1(const double &dist, double * vd1) = 0;
    virtual void urD1VD1(const double &dist, double * d1vd1) = 0;
    virtual void urD2VD1(const double &dist, double * d1vd1) = 0;

};


#endif
