#include <iostream>
#include <cmath>

#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "VMC.hpp"
#include "ConjGrad.hpp"


class QuadrExponential1D1POrbital: public WaveFunction
{
   protected:
      double _a, _b;

   public:
      QuadrExponential1D1POrbital(const double a, const double b): WaveFunction(1,1,1,2) {_a=a; _b=b;}

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
         out[0]=_a; out[1]=_b;
      }

      void samplingFunction(const double *in, double *out)
      {
         *out = -2.*(_b*(in[0]-_a)*(in[0]-_a));
      }

      double getAcceptance()
      {
         return exp(getProtoNew(0)-getProtoOld(0));
      }

      double d1(const int &i, const double *in)
      {
         return (-2.*_b*(in[0]-_a) );
      }

      double d2(const int &i, const int &j, const double *in)
      {
         return ( -2.*_b + (-2.*_b*(in[0]-_a))*(-2.*_b*(in[0]-_a)) ) ;
      }

      double vd1(const int &i, const double *in)
      {
         if (i==0)
         {
            return (2.*_b*(in[0]-_a));
         } else if (i==1)
         {
            return (-(in[0]-_a)*(in[0]-_a));
         } else
         {
            using namespace std;
            cout << "ERRORE vd1 QuadrExponential! " << endl;
            return 0.;
         }
      }
};

class Gaussian1D1POrbital: public WaveFunction
{
   protected:
      double _b;

   public:
      Gaussian1D1POrbital(const double b): WaveFunction(1,1,1,1) {_b=b;}

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

      double d1(const int &i, const double *in)
      {
         return -2.*_b*(*in);
      }

      double d2(const int &i, const int &j, const double *in)
      {
         if (_b<0.) return -1.;
         return -2.*_b+4.*_b*_b*(*in)*(*in);
      }

      double vd1(const int &i, const double *in)
      {
         return (-(*in)*(*in));
      }
};

class HarmonicOscialltor1D1P: public Hamiltonian
{
   protected:
      double _w;
   public:
      HarmonicOscialltor1D1P(const double w, WaveFunction * wf): Hamiltonian(1, 1, wf) {_w=w;}
      double localPotentialEnergy(const double *r)
      {
         return (0.5*_w*_w*(*r)*(*r));
      }
};


int main(){
   using namespace std;

   const long NMC = 10000l;
   double ** irange = new double*[1];
   *irange = new double[2];
   irange[0][0] = -25.;
   irange[0][1] = 25.;
   
   Gaussian1D1POrbital * gauss = new Gaussian1D1POrbital(1.2);
   HarmonicOscialltor1D1P * harm_osc = new HarmonicOscialltor1D1P(1.,gauss);
   QuadrExponential1D1POrbital * qexp = new QuadrExponential1D1POrbital(1.0,1.1);
   HarmonicOscialltor1D1P * harm_osc2 = new HarmonicOscialltor1D1P(1.,qexp);


   cout << endl << " - - - EVALUATION OF ENERGY - - - " << endl << endl;
   double * b1 = new double;
   gauss->getVP(b1);
   cout << "Wave Function b     = " << *b1 << endl;
   delete b1;
   VMC * vmc = new VMC(gauss,harm_osc,NMC,10l*NMC);
   vmc->getEnergyMCI()->setIRange(irange);
   double * energy = new double[4];
   double * d_energy = new double[4];
   double * gradE = new double[1];
   double * dgradE = new double[1];
   vmc->computeEnergy(NMC,energy,d_energy);
   cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;
   vmc->computeEnergyGradient(10l*NMC,gradE,dgradE);
   cout << "Energy Gradient     = " << gradE[0] <<  " +- " << dgradE[0] << endl << endl;

   cout << endl << " - - - ONE-DIMENSIONAL MINIMIZATION - - - " << endl << endl;
   double * b = new double;
   gauss->getVP(b);
   cout << "Wave Function b     = " << *b << endl;
   cout << "Conjugate Gradient Minimization ... " << endl;
   vmc->conjugateGradientOptimization();
   gauss->getVP(b);
   cout << "Wave Function b     = " << *b << endl << endl;
   vmc->computeEnergy(NMC,energy,d_energy);
   cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;
   delete b;

   cout << endl << " - - - MULTIDIMENSIONAL MINIMIZATION - - - " << endl << endl;
   double * a1 = new double[2];
   qexp->getVP(a1);
   cout << "Wave Function a     = " << a1[0] << endl;
   cout << "Wave Function b     = " << a1[1] << endl;
   delete[] a1;
   VMC * vmc2 = new VMC(qexp,harm_osc2,NMC,10l*NMC);
   vmc2->getEnergyMCI()->setIRange(irange);
   vmc2->computeEnergy(NMC,energy,d_energy);
   cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;
   delete[] gradE; gradE = new double[2];
   delete[] dgradE; dgradE = new double[2];
   vmc2->computeEnergyGradient(10l*NMC,gradE,dgradE);
   cout << "Energy Gradient 1   = " << gradE[0] <<  " +- " << dgradE[0] << endl;
   cout << "Energy Gradient 2   = " << gradE[1] <<  " +- " << dgradE[1] << endl << endl;
   double * a = new double[2];
   qexp->getVP(a);
   cout << "Wave Function a     = " << a[0] << endl;
   cout << "Wave Function b     = " << a[1] << endl;
   cout << "Conjugate Gradient Minimization ... " << endl;
   vmc2->conjugateGradientOptimization();
   qexp->getVP(a);
   cout << "Wave Function a     = " << a[0] << endl;
   cout << "Wave Function b     = " << a[1] << endl;
   vmc2->computeEnergy(NMC,energy,d_energy);
   cout << "Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;
   delete[] a;

   delete[] gradE;
   delete[] dgradE;
   delete[] energy;
   delete[] d_energy;

   delete vmc;
   delete vmc2;
   delete harm_osc;
   delete harm_osc2;
   delete gauss;
   delete qexp;

   delete[] *irange;
   delete[] irange;

   return 0;
}


