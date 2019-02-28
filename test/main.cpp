#include <cmath>
#include <iostream>

#include "vmc/Hamiltonian.hpp"
#include "vmc/MPIVMC.hpp"
#include "vmc/VMC.hpp"
#include "vmc/WaveFunction.hpp"

#include "TestVMCFunctions.hpp"


int main(){
    using namespace std;

    MPIVMC::Init();

    const int NMC = 4000l;

    auto * gauss = new Gaussian1D1POrbital(1.2);
    auto * harm_osc = new HarmonicOscillator1D1P(1., gauss);
    // QuadrExponential1D1POrbital * qexp = new QuadrExponential1D1POrbital(1.0, 1.1);
    // HarmonicOscillator1D1P * harm_osc2 = new HarmonicOscillator1D1P(1., qexp);

    cout << endl << " - - - EVALUATION OF ENERGY - - - " << endl << endl;
    double b1;
    gauss->getVP(&b1);
    cout << "Wave Function b     = " << b1 << endl;
    VMC * vmc = new VMC(gauss, harm_osc);
    vmc->getMCI()->setIRange(-10., 10.);
    double energy[4];
    double d_energy[4];
    vmc->computeVariationalEnergy(NMC, energy, d_energy, false, false);
    cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;

    // cout << endl << " - - - ONE-DIMENSIONAL MINIMIZATION - - - " << endl << endl;
    // double * b = new double;
    // gauss->getVP(b);
    // cout << "Wave Function b     = " << *b << endl;
    // cout << "Conjugate Gradient Minimization ... " << endl;
    // vmc->conjugateGradientOptimization(NMC, 4*NMC);
    // gauss->getVP(b);
    // cout << "Wave Function b     = " << *b << endl << endl;
    // vmc->computeVariationalEnergy(NMC, energy, d_energy);
    // cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    // cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    // cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    // cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;
    // delete b;

    //cout << endl << " - - - MULTIDIMENSIONAL MINIMIZATION - - - " << endl << endl;
    //double * a1 = new double[2];
    //qexp->getVP(a1);
    //cout << "Wave Function a     = " << a1[0] << endl;
    //cout << "Wave Function b     = " << a1[1] << endl;
    //delete[] a1;
    //VMC * vmc2 = new VMC(qexp, harm_osc2);
    //vmc2->computeVariationalEnergy(NMC, energy, d_energy);
    //cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    //cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    //cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    //cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;
    //double * a = new double[2];
    //qexp->getVP(a);
    //cout << "Wave Function a     = " << a[0] << endl;
    //cout << "Wave Function b     = " << a[1] << endl;
    //cout << "Conjugate Gradient Minimization ... " << endl;
    //vmc2->conjugateGradientOptimization(NMC, 4*NMC);
    //qexp->getVP(a);
    //cout << "Wave Function a     = " << a[0] << endl;
    //cout << "Wave Function b     = " << a[1] << endl;
    //vmc2->computeVariationalEnergy(NMC, energy, d_energy);
    //cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
    //cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
    //cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
    //cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;
    //delete[] a;
    //delete vmc2;

    delete vmc;

    delete harm_osc;
    // delete harm_osc2;
    delete gauss;
    // delete qexp;

    MPIVMC::Finalize();

    return 0;
}
