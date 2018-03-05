#ifndef TWO_BODY_PSEUDO_POTENTIAL
#define  TWO_BODY_PSEUDO_POTENTIAL


class TwoBodyPseudoPotential{

private:
    ParticlesDistance * _dist;
    int _npart;

public:
    TwoBodyPseudoPotential(ParticlesDistance * _dist, const int &npart){
        _nspacedim = nspacedim;
        _npart = npart;
    }
    virtual ~TwoBodyPseudoPotential();

    int getNSpaceDim(){return _dist->getNSpaceDim();}
    int getNPart(){return _npart;}


    virtual int getNVP() = 0;
    virtual void setVP(const double *vp) = 0;
    virtual void getVP(double *vp) = 0;


    virtual double u(const double & dist) = 0;
    double u(const double * r1, const double * r2){return u(_dist());}

    virtual double d1(const double & dist);
    void d1(const double * r1, const double * r2, double * grad);

    virtual double d2() = 0;
    virtual double vd1() = 0;
    virtual double d1vd1() = 0;
    virtual double d2vd1() = 0;

};


#endif
