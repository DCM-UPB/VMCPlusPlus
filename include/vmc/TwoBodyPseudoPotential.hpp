#ifndef VMC_TWOBODYPSEUDOPOTENTIAL_HPP
#define VMC_TWOBODYPSEUDOPOTENTIAL_HPP

#include "vmc/Metric.hpp"

namespace vmc
{

class TwoBodyPseudoPotential
{
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
    TwoBodyPseudoPotential(Metric * metric, const int &nvp, bool flag_vd1=false, bool flag_d1vd1=false, bool flag_d2vd1=false);
    virtual ~TwoBodyPseudoPotential();

    int getNSpaceDim(){return _metric->getNSpaceDim();}
    int getNVP(){return _nvp;}

    // -- Pair pseudopotential
    double u(const double * r1, const double * r2);

    // --- Derivatives
    // compute all the derivatives
    void computeAllDerivatives(const double * r1, const double * r2);
    // get the derivatives
    double getD1(const int &id1){return _d1[id1];}
    double getD2(const int &id2){return _d2[id2];}
    double getVD1(const int &ivd1){return _vd1[ivd1];}
    double getD1VD1(const int &id1, const int &ivd1){return _d1vd1[id1][ivd1];}
    double getD2VD1(const int &id2, const int &ivd1){return _d2vd1[id2][ivd1];}
    // has the variational derivatives?
    bool hasVD1(){return _flag_vd1;}
    bool hasD1VD1(){return _flag_d1vd1;}
    bool hasD2VD1(){return _flag_d2vd1;}


    // --- Methods that must be implemented
    // manage variational parameters
    virtual void setVP(const double *vp) = 0;
    virtual void getVP(double *vp) = 0;
    // functions of the distance that define the pseudopotential
    virtual double ur(const double &r) = 0;                           // e.g. -b/r^5
    virtual double urD1(const double &r) = 0;                         // e.g. 5*b/r^6
    virtual double urD2(const double &r) = 0;                         // e.g. -30*b/r^7
    virtual void urVD1(const double &r, double * vd1) = 0;            // e.g. -1/r^5
    virtual void urD1VD1(const double &r, double * d1vd1) = 0;        // e.g. 5/r^6
    virtual void urD2VD1(const double &r, double * d1vd1) = 0;        // e.g. -30/r^7
};
} // namespace vmc

#endif
