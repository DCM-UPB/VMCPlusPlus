#include "ParticleArrayManager.hpp"

#include <assert.h>


int main(){
    const int NPART = 3;
    const int NSPACEDIM = 3;
    const int RANDOMINT = 37;

    ParticleArrayManager * pam = new ParticleArrayManager(NSPACEDIM);


    // create the initial multi-particle array
    double * x = new double[NPART*NSPACEDIM];
    for (int i=0; i<NPART; ++i){
        for (int j=0; j<NSPACEDIM; ++j){
            x[i*NPART + j] = RANDOMINT*i + j;
        }
    }



    // check getParticleArray
    double * vec;
    for (int i=0; i<NPART; ++i){
        vec = pam->getParticleArray(x, i);
        for (int j=0; j<NSPACEDIM; ++j){
            assert( vec[j] == x[i*NPART + j] );
            assert( x[i*NPART + j] == RANDOMINT*i + j );    // array x was not modified
        }
    }



    // check setParticleArray
    double * newx = new double[NSPACEDIM];
    for (int i=0; i<NSPACEDIM; ++i) newx[i] = 2*RANDOMINT*i;
    double * oldx = new double[NSPACEDIM];
    for (int i=0; i<NPART; ++i){
        // store the original array values
        vec = pam->getParticleArray(x, i);
        for (int j=0; j<NSPACEDIM; ++j) oldx[j] = vec[j];
        // change the values for the particle i
        pam->setParticleArray(x, i, newx);
        // check that the values in the array x are as expected
        for (int i2=0; i2<NPART; ++i2){
            if (i2 == i){
                // new values
                for (int j=0; j<NSPACEDIM; ++j){
                    assert( vec[j] == newx[j]);
                }
            } else {
                // old values
                for (int j=0; j<NSPACEDIM; ++j){
                    assert( x[i2*NPART + j] == RANDOMINT*i2 + j );
                }
            }
        }
        // set x to its original values
        pam->setParticleArray(x, i, oldx);
        // check that x is set back to the original values
        for (int i2=0; i2<NPART; ++i2){
            for (int j=0; j<NSPACEDIM; ++j){
                assert( x[i2*NPART + j] == RANDOMINT*i2 + j );
            }
        }
    }



    // check addArrayToParticleArray
    for (int i=0; i<NPART; ++i){
        // store the original array values
        vec = pam->getParticleArray(x, i);
        for (int j=0; j<NSPACEDIM; ++j) oldx[j] = vec[j];
        // change the values for the particle i
        pam->addArrayToParticleArray(x, i, newx);
        // check that the values in the array x are as expected
        for (int i2=0; i2<NPART; ++i2){
            if (i2 == i){
                // new values
                for (int j=0; j<NSPACEDIM; ++j){
                    assert( vec[j] == oldx[j] + newx[j]);
                }
            } else {
                // old values
                for (int j=0; j<NSPACEDIM; ++j){
                    assert( x[i2*NPART + j] == RANDOMINT*i2 + j );
                }
            }
        }
        // set x to its original values
        for (int j=0; j<NSPACEDIM; ++j) newx[j] = - newx[j];
        pam->addArrayToParticleArray(x, i, newx);
        for (int j=0; j<NSPACEDIM; ++j) newx[j] = - newx[j];
        // check that x is set back to the original values
        for (int i2=0; i2<NPART; ++i2){
            for (int j=0; j<NSPACEDIM; ++j){
                assert( x[i2*NPART + j] == RANDOMINT*i2 + j );
            }
        }
    }







    delete[] oldx;
    delete[] newx;
    delete[] x;
    delete pam;

    return 0;
}
