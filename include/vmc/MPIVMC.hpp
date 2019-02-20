#ifndef MPIVMC
#define MPIVMC

#include "mci/MPIMCI.hpp"
#include "mci/MCIntegrator.hpp"

#include <iostream>
#include <string>

namespace MPIVMC
{
    int MyRank()
    {
        #if USE_MPI==1
        return MPIMCI::myrank();
        #else
        return 0;
        #endif
    }

    int Size()
    {
        #if USE_MPI==1
        return MPIMCI::size();
        #else
        return 1;
        #endif
    }

    int Init()
    {
        #if USE_MPI==1
        return MPIMCI::init();
        #else
        return 0;
        #endif
    }

    void SetSeed(MCI * const mci, const std::string &filename, const int &offset = 0)
    {
        #if USE_MPI==1
        MPIMCI::setSeed(mci, filename, offset);
        #else
        std::cout << "Warning: In non-MPI mode MPIVMC::SetSeed ignores the seedfile and sets the seed to offset (default 0)." << std::endl;
        mci->setSeed(offset);
        #endif
    }

    void Integrate(MCI * const mci, const long &Nmc, double * average, double * error, const bool findMRT2step=true, const bool initialdecorrelation=true, const bool randomizeWalkers = false)
    {
        if (randomizeWalkers) mci->newRandomX();
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
