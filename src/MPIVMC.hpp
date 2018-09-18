#ifndef MPIVMC
#define MPIVMC

#include "MCIntegrator.hpp"

#if USE_MPI==1
#include "MPIMCI.hpp"
#else
namespace MPIMCI // create replacement for actual MPIMCI namespace
{
    void integrate(MCI * const mci, const long &Nmc, double * average, double * error, bool findMRT2step=true, bool initialdecorrelation=true)
    {
        mci->integrate(Nmc, average, error, findMRT2step, initialdecorrelation);
    }
}
#endif
#endif
