#ifndef MPIVMC
#define MPIVMC

#if USE_MPI==1
#include "mci/MPIMCI.hpp"
#else
#include "mci/MCIntegrator.hpp"
#endif

#include <iostream>
#include <string>

namespace MPIVMC
{
    int Rank()
    {
        #if USE_MPI==1
        return MPIMCI::rank();
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

    void Integrate(MCI * const mci, const long &Nmc, double * average, double * error, int NfindMRT2stepIterations, int NdecorrelationSteps, bool randomizeWalkers = false, int nblocks=-1)
    {
        if (randomizeWalkers) {
            mci->newRandomX();
        }
        #if USE_MPI==1
        MPIMCI::integrate(mci, Nmc, average, error, NfindMRT2stepIterations, NdecorrelationSteps, nblocks < 0 ? 16 : nblocks); // if compiling with USE_MPI, set default block count to 16
        #else
        mci->integrate(Nmc, average, error, NfindMRT2stepIterations, NdecorrelationSteps, nblocks < 0 ? 0 : nblocks);
        #endif
    }

    void Integrate(MCI * const mci, const long &Nmc, double * average, double * error, bool findMRT2step=true, bool initialdecorrelation=true, int nblocks=-1)
    {
        #if USE_MPI==1
        // if compiling with USE_MPI, set fixed defaults (totaling 5000 quick no-sample steps)
        int stepsMRT2 = findMRT2step ? 25 : 0;
        int stepsDecorr = initialdecorrelation ? 2500 : 0;
        Integrate(mci, Nmc, average, error, stepsMRT2, stepsDecorr, true, nblocks < 0 ? 16 : nblocks); // if compiling with USE_MPI, also use random initial walker positions and fixed blocking
        #else
        int stepsMRT2 = findMRT2step ? -1 : 0;
        int stepsDecorr = initialdecorrelation ? -1 : 0;
        Integrate(mci, Nmc, average, error, stepsMRT2, stepsDecorr, false, nblocks < 0 ? 0 : nblocks);
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
