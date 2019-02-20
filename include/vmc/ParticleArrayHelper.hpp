#ifndef PARTICLE_ARRAY_HELPER
#define PARTICLE_ARRAY_HELPER

/*
Tools for managing 1D arrays which represent the matrix:
    x[i][j] = coordinate j of the particle i

It assumes that the array x is built in this way:
        x[ i*nspacedim + j ]

Example: npart=2 nspacedim=3   ->   1st particle = {x1, x2, x3}   2nd particle = {y1, y2, y3}
    x = {x1, x2, x3, y1, y2, y3}

IMPORTANT: the out-of boundary check (for the index i, which might exceed npart) is not made for perfomance reasons.
Since this class most likely will be used within the library, unit tests on the classes which make uses of this class should rule out errors.
*/


class ParticleArrayHelper
{
private:
    int _nspacedim;

public:
    explicit ParticleArrayHelper(const int &nspacedim){
        _nspacedim = nspacedim;
    }

    double * getParticleArray(double * x, const int &i);
    const double * getParticleArray(const double * x, const int &i);

    void setParticleArray(double * x, const int &i, const double * newx);

    void addArrayToParticleArray(double * x, const int &i, const double * toadd);
};


#endif
