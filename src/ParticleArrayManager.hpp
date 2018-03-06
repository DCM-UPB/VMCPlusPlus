#ifndef PARTICLE_ARRAY_MANAGER
#define PARTICLE_ARRAY_MANAGER

/*
Assumes that the array x is built in this way:
        x[ i*npart + j ]
i.e. i is the particle index and j the space index.

Example: npart=2 nspacedim=3   ->   1st particle = {x1, x2, x3}   2nd particle = {y1, y2, y3}
    x = {x1, x2, x3, y1, y2, y3}

TODO

IMPORTANT: the out-of boundary checks (for the index i, which might exceed npart) is not made for perfomance reasons.
Since this class most likely will be used within the library, unit tests on the classes that makes uses of this class should rule out errors.
*/


class ParticleArrayManager{

private:
    int _nspacedim;

public:
    ParticleArrayManager(const int &nspacedim){
        _nspacedim = nspacedim;
    }

    double * getParticleArray(double * x, const int &i);
    const double * getParticleArray(const double * x, const int &i);

    void setParticleArray(double * x, const int &i, const double * newx);

    void addArrayToParticleArray(double * x, const int &i, const double * toadd);

};


#endif
