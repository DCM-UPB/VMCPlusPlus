#ifndef PARTICLE_DISTANCE
#define PARTICLE_DISTANCE


class ParticlesDistance{

private:
    int _nspacedim;

public:
    ParticlesDistance(const int &nspacedim){
        _nspacedim = nspacedim;
    }


    int getNSpaceDim(){return _nspacedim;}



    virtual double dist(const double * r1, const double * r2) = 0;

    virtual void distD1(const double * r1, const double * r2, double * out) = 0;

    virtual void distD2(const double * r1, const double * r2, double * out) = 0;

};


#endif
