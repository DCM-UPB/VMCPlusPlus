#ifndef VMC_MPIVMC_HPP
#define VMC_MPIVMC_HPP

#include "mci/MCIntegrator.hpp"
#include "mci/MPIMCI.hpp"

#include <iostream>
#include <string>

namespace MPIVMC
{
inline int MyRank()
{
#if USE_MPI == 1
    return MPIMCI::myrank();
#else
    return 0;
#endif
}

inline int Size()
{
#if USE_MPI == 1
    return MPIMCI::size();
#else
    return 1;
#endif
}

inline int Init()
{
#if USE_MPI == 1
    return MPIMCI::init();
#else
    return 0;
#endif
}

inline void SetSeed(mci::MCI &mci, const std::string &filename, int offset = 0)
{
#if USE_MPI == 1
    MPIMCI::setSeed(mci, filename, offset);
#else
    std::cout << "Warning: In non-MPI mode MPIVMC::SetSeed ignores the seedfile and sets the seed to offset (default 0)." << std::endl;
    mci.setSeed(offset);
#endif
}

inline void Integrate(mci::MCI &mci, int Nmc, double * average, double * error, bool findMRT2step = true, bool initialdecorrelation = true, bool randomizeWalkers = false)
{
    if (randomizeWalkers) { mci.newRandomX(); }
#if USE_MPI == 1
    MPIMCI::integrate(mci, Nmc, average, error, findMRT2step, initialdecorrelation);
#else
    mci.integrate(Nmc, average, error, findMRT2step, initialdecorrelation);
#endif
}

inline void Finalize()
{
#if USE_MPI == 1
    MPIMCI::finalize();
#endif
}
}  // namespace MPIVMC

#endif
