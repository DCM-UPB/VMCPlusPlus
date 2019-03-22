#include "vmc/SymmetrizerWaveFunction.hpp"
#include "vmc/PairSymmetrizerWaveFunction.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>

#include "TestVMCFunctions.hpp"


int main(){
    using namespace std;

    // constants
    const int NSPACEDIM = 1;
    const int NPART = 3;
    const int NTOTALDIM = NPART*NSPACEDIM; 
    const double GAUSS_EXPF = 0.9;
    const double DX = 0.0001;
    const double TINY = 0.02;
    const double SUPERTINY = 0.0000001;

    // random generator
    random_device rdev;
    mt19937_64 rgen;
    uniform_real_distribution<double> rd;
    rgen = mt19937_64(rdev());
    rgen.seed(18984687);
    rd = uniform_real_distribution<double>(-0.05, 0.05);

    // define a non-symmetric wavefunction
    const double ai_nosym[NPART] = {0.5, -0.25, 0.0};
    auto * phi_nosym = new QuadrExponential1DNPOrbital(NPART, ai_nosym, GAUSS_EXPF);

    // fully (anti-)symmetrized wfs
    auto * phi_sym = new SymmetrizerWaveFunction(phi_nosym, false);
    auto * phi_asym = new SymmetrizerWaveFunction(phi_nosym, true); // antisymmetric version

    // particles position and all permutations
    double * xp[6];
    for (auto & x : xp) { x = new double[NTOTALDIM]; }

    xp[0][0] = 0.2; xp[0][1] = -0.5; xp[0][2] = 0.7;
    xp[1][0] = xp[0][0]; xp[1][1] = xp[0][2]; xp[1][2] = xp[0][1];
    xp[2][0] = xp[0][1]; xp[2][1] = xp[0][0]; xp[2][2] = xp[0][2];
    xp[3][0] = xp[0][1]; xp[3][1] = xp[0][2]; xp[3][2] = xp[0][0];
    xp[4][0] = xp[0][2]; xp[4][1] = xp[0][0]; xp[4][2] = xp[0][1];
    xp[5][0] = xp[0][2]; xp[5][1] = xp[0][1]; xp[5][2] = xp[0][0];

    bool isOdd[6] = {false, true, true, false, false, true};

    // compute the new protovalues
    double protov_nosym[2], protov_sym[2], protov_asym[2];
    
    for (int i=0; i<6; ++i) {
        size_t offset = (i==0 ? 0 : 1);

        phi_nosym->protoFunction(xp[i], protov_nosym + offset);
        phi_sym->protoFunction(xp[i], protov_sym + offset);
        phi_asym->protoFunction(xp[i], protov_asym + offset);

        // cout << "Perm " << i << ": ";
        // cout << " phi_nosym " << phi_nosym->computeWFValue(protov_nosym + offset);
        // cout << " phi_sym " << phi_sym->computeWFValue(protov_sym + offset);
        // cout << " phi_asym " << phi_asym->computeWFValue(protov_asym + offset) << endl;

        if (i>0) {
            assert(protov_nosym[0] != protov_nosym[1]);
            if (isOdd[i]) { assert(fabs(protov_asym[0] + protov_asym[1]) < SUPERTINY);
            } else { assert(fabs(protov_asym[0] - protov_asym[1]) < SUPERTINY); }

            assert(fabs(protov_sym[0] - protov_sym[1]) < SUPERTINY);
        }
    }

    std::vector<WaveFunction *> wfs;
    wfs.push_back(phi_nosym);
    wfs.push_back(phi_sym);
    wfs.push_back(phi_asym);

    vector<string> names {"phi_nosym", "phi_sym", "phi_asym"};

    double x[NTOTALDIM];
    for (int i=0; i<NTOTALDIM; ++i) { x[i] = xp[0][i]; }

    double vp[1];
    vp[0] = GAUSS_EXPF;
    int cont = 0;
    for (WaveFunction * wf : wfs) {
        // cout << "Checking " << names[cont] << " ..." << endl;
        
        // verify variational parameter
        double wfvp;
        wf->getVP(&wfvp);
        assert(vp[0] == wfvp);

        // pre-compute all the derivatives analytically
        wf->computeAllDerivatives(xp[0]);

        // initial wave function
        double f, fdx, fmdx, fdvp, fdxdvp, fmdxdvp;
        double samp;
        wf->protoFunction(x, &samp); f = wf->computeWFValue(&samp);


        // --- check the first derivatives
        for (int i=0; i<NTOTALDIM; ++i){
            const double origx = x[i];
            x[i] += DX;
            wf->protoFunction(x, &samp); fdx = wf->computeWFValue(&samp);
            const double numderiv = (fdx-f)/(DX*f);

            // cout << "getD1DivByWF(" << i <<") = " << wf->getD1DivByWF(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs( (wf->getD1DivByWF(i) - numderiv)/numderiv) < TINY );

            x[i] = origx;
        }


        // --- check the second derivatives
        for (int i=0; i<NTOTALDIM; ++i){
            const double origx = x[i];
            x[i] = origx + DX;
            wf->protoFunction(x, &samp); fdx = wf->computeWFValue(&samp);
            x[i] = origx - DX;
            wf->protoFunction(x, &samp); fmdx = wf->computeWFValue(&samp);
            const double numderiv = (fdx - 2.*f + fmdx) / (DX*DX*f);

            // cout << "getD2DivByWF(" << i << ") = " << wf->getD2DivByWF(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs( (wf->getD2DivByWF(i) - numderiv)/numderiv) < TINY );

            x[i] = origx;
        }


        // -- check the first variational derivative
        for (int i=0; i<wf->getNVP(); ++i){
            const double origvp = vp[i];
            vp[i] += DX;
            wf->setVP(vp);
            wf->protoFunction(x, &samp); fdvp = wf->computeWFValue(&samp);
            const double numderiv = (fdvp - f)/(DX*f);

            // cout << "getVD1DivByWF(" << i << ") = " << wf->getVD1DivByWF(i) << endl;
            // cout << "numderiv = " << numderiv << endl << endl;
            assert( fabs( (wf->getVD1DivByWF(i) - numderiv)/numderiv ) < TINY );

            vp[i] = origvp;
            wf->setVP(vp);
        }


        // --- check the first cross derivative
        for (int i=0; i<NTOTALDIM; ++i){
            for (int j=0; j<wf->getNVP(); ++j){
                const double origx = x[i];
                const double origvp = vp[j];

                x[i] += DX;
                wf->protoFunction(x, &samp); fdx = wf->computeWFValue(&samp);

                x[i] = origx;
                vp[j] += DX;
                wf->setVP(vp);
                wf->protoFunction(x, &samp); fdvp = wf->computeWFValue(&samp);

                x[i] += DX;
                wf->protoFunction(x, &samp); fdxdvp = wf->computeWFValue(&samp);

                const double numderiv = (fdxdvp - fdx - fdvp + f)/(DX*DX*f);

                // cout << "getD1VD1DivByWF(" << i << ", " << j << ") = " << wf->getD1VD1DivByWF(i, j) << endl;
                // cout << "numderiv = " << numderiv << endl << endl;
                assert( fabs( (wf->getD1VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );

                x[i] = origx;
                vp[j] = origvp;
                wf->setVP(vp);
            }
        }


        // --- check the second cross derivative
        for (int i=0; i<NTOTALDIM; ++i){
            for (int j=0; j<wf->getNVP(); ++j){
                const double origx = x[i];
                const double origvp = vp[j];

                x[i] = origx + DX;
                wf->protoFunction(x, &samp); fdx = wf->computeWFValue(&samp);

                vp[j] = origvp + DX;
                wf->setVP(vp);
                wf->protoFunction(x, &samp); fdxdvp = wf->computeWFValue(&samp);

                x[i] = origx;
                wf->protoFunction(x, &samp); fdvp = wf->computeWFValue(&samp);

                x[i] = origx - DX;
                wf->protoFunction(x, &samp); fmdxdvp = wf->computeWFValue(&samp);

                vp[j] = origvp;
                wf->setVP(vp);
                wf->protoFunction(x, &samp); fmdx = wf->computeWFValue(&samp);

                const double numderiv = (fdxdvp - 2.*fdvp + fmdxdvp - fdx + 2.*f - fmdx)/(DX*DX*DX*f);

                // cout << "getD2VD1DivByWF(" << i << ", " << j << ") = " << wf->getD2VD1DivByWF(i, j) << endl;
                // cout << "numderiv = " << numderiv << endl << endl;
                assert( fabs( (wf->getD2VD1DivByWF(i, j)-numderiv)/numderiv ) < TINY );

                x[i] = origx;
                vp[j] = origvp;
                wf->setVP(vp);
            }
        }
        ++cont;
    }

    for (auto & x : xp) { delete [] x; }
    
    delete phi_asym;
    delete phi_sym;

    delete phi_nosym;

    return 0;
}
