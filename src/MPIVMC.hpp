#ifndef MPIVMC
#define MPIVMC

#if USE_MPI==1
#include "MPIMCI.hpp"
#else
#include "MCIntegrator.hpp"
#endif

namespace MPIVMC
{
    int Init()
    {
        #if USE_MPI==1
        return MPIMCI::init();
        #else
        return 0;
        #endif
    }

    void Integrate(MCI * const mci, const long &Nmc, double * average, double * error, bool findMRT2step=true, bool initialdecorrelation=true)
    {
        #if USE_MPI==1
        MPIMCI::integrate(mci, Nmc, average, error, findMRT2step, initialdecorrelation);
        #else
        mci->integrate(Nmc, average, error, findMRT2step, initialdecorrelation);
        #endif
    }

    void Finalize()
    {
        #if USE_MPI==1
        MPIMCI::finalize();
        #endif
    }
};

#endif
