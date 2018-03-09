#include <iostream>
#include <cmath>

#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "VMC.hpp"




class QuadrExponential1D1POrbital: public WaveFunction
{
protected:
    double _a, _b;
    double _wf_exp, _d1, _d2, _vd1_a, _vd1_b;

public:
    QuadrExponential1D1POrbital(const double a, const double b): WaveFunction(1, 1, 1, 2, true, false, false) {_a=a; _b=b;}

    void setVP(const double *in)
    {
        _a=in[0];
        //if (_a<0.01) _a=0.01;
        _b=in[1];
        //if (_b<0.01) _b=0.01;
        using namespace std;
        //cout << "change a and b! " << _a << "   " << _b << endl;
    }
    void getVP(double *out)
    {
        out[0]=_a;
        out[1]=_b;
    }

    void samplingFunction(const double *in, double *out)
    {
        *out = -2.*(_b*(in[0]-_a)*(in[0]-_a));
    }

    double getAcceptance()
    {
        return exp(getProtoNew(0)-getProtoOld(0));
    }

    void computeAllDerivatives(const double *in){
        _setD1DivByWF(0, -2.*_b*(in[0]-_a));
        _setD2DivByWF(0, -2.*_b + (-2.*_b*(in[0]-_a))*(-2.*_b*(in[0]-_a)));
        if (hasVD1()){
            _setVD1DivByWF(0, 2.*_b*(in[0]-_a));
            _setVD1DivByWF(1, -(in[0]-_a)*(in[0]-_a));
        }
    }
};



class Gaussian1D1POrbital: public WaveFunction
{
protected:
    double _b;

public:
    Gaussian1D1POrbital(const double b):
    WaveFunction(1, 1, 1, 1, false, false, false){
        _b=b;
    }

    void setVP(const double *in)
    {
        _b=*in;
        //if (_b<0.01) _b=0.01;
        using namespace std;
        //cout << "change b! " << _b << endl;
    }
    void getVP(double *out)
    {
        *out=_b;
    }

    void samplingFunction(const double *in, double *out)
    {
        *out=-2.*_b*(*in)*(*in);
    }

    double getAcceptance()
    {
        return exp(getProtoNew(0)-getProtoOld(0));
    }

    void computeAllDerivatives(const double *in){
        _setD1DivByWF(0, -2.*_b*(*in));
        _setD2DivByWF(0, -2.*_b+4.*_b*_b*(*in)*(*in));
        if (hasVD1()){
            _setVD1DivByWF(0, (-(*in)*(*in)));
        }
    }
};



class HarmonicOscillator1D1P: public Hamiltonian
{
protected:
    double _w;
public:
    HarmonicOscillator1D1P(const double w, WaveFunction * wf): Hamiltonian(1, 1, wf) {_w=w;}
    double localPotentialEnergy(const double *r)
    {
        return (0.5*_w*_w*(*r)*(*r));
    }
};


int main(){
    using namespace std;

    const long NMC = 4000l;
    double ** irange = new double*[1];
    *irange = new double[2];
    irange[0][0] = -25.;
    irange[0][1] = 25.;

    Gaussian1D1POrbital * gauss = new Gaussian1D1POrbital(1.2);
    HarmonicOscillator1D1P * harm_osc = new HarmonicOscillator1D1P(1., gauss);
    // QuadrExponential1D1POrbital * qexp = new QuadrExponential1D1POrbital(1.0, 1.1);
    // HarmonicOscillator1D1P * harm_osc2 = new HarmonicOscillator1D1P(1., qexp);

    cout << endl << " - - - EVALUATION OF ENERGY - - - " << endl << endl;
    double * b1 = new double;
    gauss->getVP(b1);
    cout << "Wave Function b     = " << *b1 << endl;
    VMC * vmc = new VMC(gauss, harm_osc);
    vmc->getMCI()->setIRange(irange);
    double * energy = new double[4];
    double * d_energy = new double[4];
    vmc->computeVariationalEnergy(NMC, energy, d_energy);
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

    delete[] energy;
    delete[] d_energy;

    delete vmc;
    delete b1;

    delete harm_osc;
    // delete harm_osc2;
    delete gauss;
    // delete qexp;

    delete[] *irange;
    delete[] irange;

    return 0;
}
