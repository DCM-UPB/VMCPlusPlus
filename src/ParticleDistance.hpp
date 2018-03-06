#ifndef PARTICLE_DISTANCE
#define PARTICLE_DISTANCE


class ParticleDistance{

private:
    int _nspacedim;

public:
    ParticleDistance(const int &nspacedim){
        _nspacedim = nspacedim;
    }
    virtual ~ParticleDistance(){}


    int getNSpaceDim(){return _nspacedim;}



    virtual double dist(const double * r1, const double * r2) = 0;

    virtual void distD1(const double * r1, const double * r2, double * out) = 0;

    virtual void distD2(const double * r1, const double * r2, double * out) = 0;

};


#endif
