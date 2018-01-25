//#include <iostream>
//#include <cmath>
//#include <math.h>
//#include <fstream>
//#include <random>
//
//#include "NoisyFunction.hpp"
//#include "ConjGrad.hpp"
//
//
//
//class Noiseless2DParabola: public NoisyFunctionWithGradient{
//public:
//   Noiseless2DParabola(): NoisyFunctionWithGradient(2){}
//   
//   void f(const double * in, double &f, double &df){
//      f = pow(in[0]-1., 2) + pow(in[1]+2., 2);   // minimum in (1, -2)
//      df = 0.;
//   }
//   
//   void grad(const double * in, double * g, double * dg){
//      g[0] = 2. * (in[0] - 1.);
//      g[1] = 2. * (in[1] + 2.);
//      dg[0] = 0.;
//      dg[1] = 0.;
//   }
//};
//
//
//
//class Noisy2DParabola: public NoisyFunctionWithGradient{
//private:
//   const double _sigma = 0.5;
//   std::random_device _rdev;
//   std::mt19937_64 _rgen;
//   std::uniform_real_distribution<double> _rd;  //after initialization (done in the constructor) can be used with _rd(_rgen)
//   
//public:
//   Noisy2DParabola(): NoisyFunctionWithGradient(2){
//      // initialize random generator
//      _rgen = std::mt19937_64(_rdev());
//      _rd = std::uniform_real_distribution<double>(-_sigma, _sigma);
//   }
//   
//   void f(const double * in, double &f, double &df){
//      f = pow(in[0]-1., 2) + pow(in[1]+2., 2);   // minimum in (1, -2)
//      df = _sigma;
//      f += _rd(_rgen);
//   }
//   
//   void grad(const double * in, double * g, double * dg){
//      g[0] = 2. * (in[0] - 1.);
//      g[1] = 2. * (in[1] + 2.);
//      dg[0] = 2.*_sigma;
//      dg[1] = 2.*_sigma;
//      g[0] += 2.*_rd(_rgen);
//      g[1] += 2.*_rd(_rgen);
//   }
//};
//
//
//
//
//int main() {
//    using namespace std;
//    
//    cout << "We want to minimize the 2D function" << endl;
//    cout << "    (x-1)^2 + (y+2)^2" << endl;
//    cout << "whose min is in (1, -2)." << endl << endl << endl;
//    
//    
//    
//    
//    cout << "we first minimize it, supposing to have no noise at all" << endl;
//    
//    Noiseless2DParabola * nlp = new Noiseless2DParabola();
//    
//    ConjGrad * cg = new ConjGrad(nlp);
//    
//    double * initpos = new double[2];
//    initpos[0] = -1.;
//    initpos[1] = -1.;
//    cg->setX(initpos);
//        
//    cg->findMin();
//    
//    cout << "The found minimum is: ";
//    cout << cg->getX(0) << "    " << cg->getX(1) << endl << endl << endl;
//    
//    
//    
//    
//    cout << "Now we repeat the minimisation adding a noise to the function and its gradient." << endl;
//    
//    Noisy2DParabola * np = new Noisy2DParabola();
//    
//    delete cg;
//    cg = new ConjGrad(np);
//    cg->setX(initpos);
//        
//    cg->findMin();
//    
//    cout << "The found minimum is: ";
//    cout << cg->getX(0) << "    " << cg->getX(1) << endl << endl;
//    
//    
//    delete np;
//    delete[] initpos;
//    delete cg;
//    delete nlp;
//    
//
//    // end
//    return 0;
//}






#include <iostream>
#include <cmath>

#include "WaveFunction.hpp"
#include "Hamiltonian.hpp"
#include "VMC.hpp"
#include "ConjGrad.hpp"



/*
Hamiltonian describing a 1-particle harmonic oscillator:
   H  =  p^2 / 2m  +  1/2 * w^2 * x^2 
*/
class HarmonicOscillator1D1P: public Hamiltonian{

protected:
   double _w;

public:
   HarmonicOscillator1D1P(const double w, WaveFunction * wf): 
      Hamiltonian(1 /*num space dimensions*/, 1 /*num particles*/, wf) {_w=w;}
   
   // potential energy
   double localPotentialEnergy(const double *r)
   {
      return (0.5*_w*_w*(*r)*(*r));
   }
};



/*
Trial Wave Function for 1 particle in 1 dimension, that uses two variational parameters: a and b.
   Psi  =  exp( -b * (x-a)^2 )
Notice that the corresponding probability density (sampling function) is Psi^2.
*/
class QuadrExponential1D1POrbital: public WaveFunction{
   protected:
      double _a, _b;

   public:
      QuadrExponential1D1POrbital(const double a, const double b): 
         WaveFunction(1 /*num space dimensions*/, 1 /*num particles*/, 1 /*num wf components*/, 2 /*num variational parameters*/) {_a=a; _b=b;}
      
      void setVP(const double *in){
         _a=in[0];
         _b=in[1];
      }
      
      void getVP(double *out){
         out[0]=_a; out[1]=_b;
      }

      void samplingFunction(const double *x, double *out){
         /*
         Compute the sampling function proto value, used in getAcceptance()
         */
         *out = -2.*(_b*(x[0]-_a)*(x[0]-_a));
      }

      double getAcceptance(){
         /*
         Compute the acceptance probability
         */
         return exp(getProtoNew(0)-getProtoOld(0));
      }

      double d1(const int &i, const double *x){
         /*
         Compute:    d/dx_i log(Psi(x))
         */
         return (-2.*_b*(x[0]-_a) );
      }

      double d2(const int &i, const int &j, const double *x){
         /*
         Compute:    d^2/dx_i^2 log(Psi(x))
         */
         return ( -2.*_b + (-2.*_b*(x[0]-_a))*(-2.*_b*(x[0]-_a)) ) ;
      }

      double vd1(const int &i, const double *x){
         /*
         Compute:    d/dalpha_i log(Psi(x))
         where alpha are the variational parameters, in our case an array of dimension 1: alpha = (b)
         */
         if (i==0){
            return (2.*_b*(x[0]-_a));
         } else if (i==1){
            return (-(x[0]-_a)*(x[0]-_a));
         } else{
            using namespace std;
            cout << "ERRORE vd1 QuadrExponential! " << endl;
            return 0.;
         }
      }
};




int main(){
   using namespace std;
   
   // Declare some trial wave functions
   QuadrExponential1D1POrbital * psi = new QuadrExponential1D1POrbital(-0.5, 1.0);
   
   // Declare an Hamiltonian
   // We use the harmonic oscillator with w=1 and w=2
   const double w1 = 1.;
   HarmonicOscillator1D1P * ham1 = new HarmonicOscillator1D1P(w1, psi);
   const double w2 = 2.;
   HarmonicOscillator1D1P * ham2 = new HarmonicOscillator1D1P(w2, psi);
   
   
   
   cout << endl << " - - - WAVE FUNCTION OPTIMIZATION - - - " << endl << endl;
   
   VMC * vmc; // VMC object we will resuse
   const long E_NMC = 4000l; // MC samplings to use for computing the energy
   const long G_NMC = 10000l; // MC samplings to use for computing the energy gradient
   double * energy = new double[4]; // energy
   double * d_energy = new double[4]; // energy error bar
   double * vp = new double[psi->getNVP()];
   
   
   
   // Case 1
   cout << "-> ham1:    w = " << w1 << endl << endl;
   vmc = new VMC(psi, ham1);
   
   cout << "   Initial Wave Function parameters:" << endl;
   psi->getVP(vp);
   cout << "       a = " << vp[0] << endl;
   cout << "       b = " << vp[1] << endl;
   
   cout << "   Starting energy:" << endl;
   vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
   cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;
   
   cout << "   Optimization . . ." << endl;
   vmc->conjugateGradientOptimization(E_NMC, G_NMC);
   cout << "   . . . Done!" << endl << endl;
   
   cout << "   Optimized Wave Function parameters:" << endl;
   psi->getVP(vp);
   cout << "       a = " << vp[0] << endl;
   cout << "       b = " << vp[1] << endl;
   
   cout << "   Optimized energy:" << endl;
   vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
   cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl << endl;
   
   
   
   // Case 2
   cout << "-> ham2:    w = " << w2 << endl << endl;
   delete vmc;
   vmc = new VMC(psi, ham2);
   
   cout << "   Initial Wave Function parameters:" << endl;
   psi->getVP(vp);
   cout << "       a = " << vp[0] << endl;
   cout << "       b = " << vp[1] << endl;
   
   cout << "   Starting energy:" << endl;
   vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
   cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;
   
   cout << "   Optimization . . ." << endl;
   vmc->conjugateGradientOptimization(E_NMC, G_NMC);
   cout << "   . . . Done!" << endl << endl;
   
   cout << "   Optimized Wave Function parameters:" << endl;
   psi->getVP(vp);
   cout << "       a = " << vp[0] << endl;
   cout << "       b = " << vp[1] << endl;
   
   cout << "   Optimized energy:" << endl;
   vmc->computeVariationalEnergy(E_NMC, energy, d_energy);
   cout << "       Total Energy        = " << energy[0] << " +- " << d_energy[0] << endl;
   cout << "       Potential Energy    = " << energy[1] << " +- " << d_energy[1] << endl;
   cout << "       Kinetic (PB) Energy = " << energy[2] << " +- " << d_energy[2] << endl;
   cout << "       Kinetic (JF) Energy = " << energy[3] << " +- " << d_energy[3] << endl << endl;
   

   
   
   delete[] vp;
   delete[] d_energy;
   delete[] energy;
   delete vmc;

   

   return 0;
}


