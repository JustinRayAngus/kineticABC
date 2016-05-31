/***
 * 
 * electron-energy-distribution function class
 *
***/

#ifndef EEDF_h
#define EEDF_h

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>

#include "Gas.h"
#include "json/json.h"
#include "energyGrid.h"
#include "HDF5dataFile.h"

using namespace std;

class EEDF
{
public:
   double zeroMom, Te0;   // zero momentum, Te;
   string type0;          // initial type of EEDF
   vector<double> F0, F0old, F0half; // EEDF
   vector<double> W, D, Flux;  // Energy space adv, diff, and flux at cell-edge
   vector<double> ExcS; // electronic-excitation source term at cell-center
   
   double dtStable;  // stable time-step
   
   void initialize(const energyGrid&, const Json::Value&, HDF5dataFile&);
   void computeFlux(const Gas&, const energyGrid&, const double&);
   void computeExcS(const Gas&, const energyGrid&);
   void advanceF0(const energyGrid&, const double&);
   
private:
   void setTriDiagCoeffs(vector<double>&, vector<double>&, vector<double>&, vector<double>&,
                         const energyGrid&, const double&);
   void IntExp(double&, const double&, const double&, const double&, const int&);
};

void EEDF::initialize(const energyGrid& Egrid, const Json::Value& root, HDF5dataFile& dataFile)
{
   F0.assign(Egrid.nE,0.0);
   Flux.assign(Egrid.nE+1,0.0);
   ExcS.assign(Egrid.nE,0.0);
   W.assign(Egrid.nE+1,0.0);
   D.assign(Egrid.nE+1,0.0);
   const Json::Value defValue; // used for default reference
   const Json::Value EEDF = root.get("EEDF",defValue);
   if(EEDF.isObject()) {
      printf("Initializing EEDF ...\n");
      Json::Value Te = EEDF.get("Te",defValue);
      Json::Value type = EEDF.get("type",defValue);
      if(Te == defValue || type == defValue) {
         printf("ERROR: Te or type not declared in input file\n");
         exit (EXIT_FAILURE);
      } 
      Te0 = Te.asDouble();
      if(Te0 < 0.0) {
         printf("ERROR: Te is not a positive value in input file\n");
         exit (EXIT_FAILURE);
      }
      type0 = type.asString();
      if(type0=="maxwellian" || type0=="Maxwellian") {
         cout << "Initial EEDF is Maxwellian with Te = " << Te0 << endl;
         zeroMom = 0;
         for (auto n=0; n<Egrid.nE; n++) {
            F0[n] = 2.0/sqrt(Te0*3.14159)/Te0*exp(-Egrid.Ecc[n]/Te0);
            zeroMom = zeroMom + sqrt(Egrid.Ecc[n])*F0[n]*Egrid.dE;
         }
      }
      else {
         cout << "Initial EEDF type = " << type0 << " is not valid " << endl;
         exit (EXIT_FAILURE);
      }
   }
   else {
      cout << "value for key \"Egrid\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   F0old = F0;
   F0half = F0;
   
   dataFile.add(F0, "F0", 1);  
   dataFile.add(Flux, "Flux", 1); //Flux = WF - D*dF/dE
   dataFile.add(W, "W", 1);  // advection coefficient
   dataFile.add(D, "D", 1);  // diffusion coefficient
   dataFile.add(ExcS, "ExcS", 1); // excitation source term
   dataFile.add(zeroMom, "zeroMom", 1); // should remain at one always !!! 
   dataFile.add(Te0, "Te0", 1); 
   cout << endl;  
}

void EEDF::computeFlux(const Gas& gas, const energyGrid& Egrid, const double& EVpm)
{
   const double Ng = gas.Ng;
   const double mM = gas.mM;
   const double econst = 1.6022e-19;
   const double meconst = 9.1094e-31;
   const double gamma = sqrt(2.0*econst/meconst);
   const int nE = Egrid.nE; // number of cell-center points
   
   // compute advection and diffusion coefficients
   //
   vector<double> PecNum, dtMaxW, dtMaxD, dtMaxWD;
   PecNum.assign(nE+1,0.0);     // Peclet number, |W*dE/D| < 2, 
   dtMaxW.assign(nE+1,1E20);    // Courant: coura=|W|*dt/(dE*sqrt(E)) < 1
   dtMaxD.assign(nE+1,1E20);    // Diff:    alpha=|U|*dt/(dE^2*sqrt(E)) < 1/2
   dtMaxWD.assign(nE+1,1E20);   // coura^2/alpha < 2
   
   for (auto n=1; n<nE+1; n++) { // dont need values at E=0
      W[n] = -gamma*Ng*2.0*mM*pow(Egrid.Ece[n],2.0)*gas.Qelm[n];
      D[n] = gamma*EVpm*EVpm/3.0/Ng/gas.Qmom[n]*Egrid.Ece[n];
      PecNum[n] = W[n]*Egrid.dE/D[n];
      dtMaxW[n] = Egrid.dE*sqrt(Egrid.Ece[n])/abs(W[n]);
      dtMaxD[n] = 0.5*Egrid.dE*Egrid.dE*sqrt(Egrid.Ece[n])/D[n];
      dtMaxWD[n] = 2.0*D[n]/W[n]/W[n];
   }
   
   // make sure grid is refined enough for results to be physical
   //
   double PecMax = *max_element(begin(PecNum), end(PecNum));
   double PecMin = *min_element(begin(PecNum), end(PecNum));
   if(PecMax >= 2 || abs(PecMin) >=2) {
      cout << "Maximum Peclet number = " << max(PecMax,abs(PecMin)) << " >=2 " << endl;
      cout << "Grid is not refined enough for results to be physical !!!" << endl;
      exit(EXIT_FAILURE);
   }
   
   // compute stable time step for explicit integration
   //
   double dtW = *min_element(begin(dtMaxW), end(dtMaxW));
   double dtD = *min_element(begin(dtMaxD), end(dtMaxD));
   double dtWD = *min_element(begin(dtMaxWD), end(dtMaxWD));
   //cout << "Stable time step for advection: " << dtW << endl;
   //cout << "Stable time step for diffusion: " << dtD << endl;
   //cout << "Stable time step for advection-diffusion: " << dtWD << endl;
   dtStable = min(min(dtW,dtD),dtWD);
      
   // compute Flux = W*F0 - D*dF0/dE
   // boundary conditions are Flux=0 at E=0 and E=Emax
   // E=0 BC is physical and E=Emax BC is to conserve zero moment
   // F0 and dF0/dE are finite at E=0 => W = D = 0 at E=0
   // W*F0=D*DF0/dE at E=Emax
   //
   Flux[0] = 0;
   for (auto n=1; n<nE; n++) { 
      Flux[n] = W[n]*(F0half[n]+F0half[n-1])/2.0 
              - D[n]*(F0half[n]-F0half[n-1])/Egrid.dE;
   }
   Flux[nE] = 0; 
}


void EEDF::computeExcS(const Gas& gas, const energyGrid& Egrid)
{
   const double Ng = gas.Ng;
   const vector<vector<double>> Qexc = gas.Qexc;
   const vector<double> Uexc = gas.Uexc;
   const double econst = 1.6022e-19;
   const double meconst = 9.1094e-31;
   const double gamma = sqrt(2.0*econst/meconst);
   const int nE = Egrid.nE; // number of cell-center points
   
   // loop through each reaction
   //
   double thisU, PL, PU, Esub, deltaE;
   double a, b;
   // deltaQ;
   //double c, d, g, exp1, exp2;
   int thism;
   vector<double> thisQ;
   int numReacs = Uexc.size();
   fill(ExcS.begin(),ExcS.end(),0.0); // reset to zeros;
   for (auto n=0; n<numReacs; n++) {
      thisU = Uexc[n];
      thisQ = Qexc[n];
      //cout << thisQ.size() << endl;
      thism = 1;
      for (auto j=0; j<nE; j++) {
         if( Egrid.Ece[j] < thisU && thisU <= Egrid.Ece[j+1] ) {
         
            // this cell contains thisU
            //
            a = thisU;
            b = Egrid.Ece[j+1];
            if(a==b) a = Egrid.Ece[j]; // prevent divide by zero below using method 2 and 3
 //           deltaQ = thisQ[j+1]-thisQ[j]; // method 2
            PU = Ng*gamma*F0half[j]*(thisQ[j]+thisQ[j+1])/2.0*(b*b-a*a)/2.0; // method 1
 //         PU = Ng*gamma*F0half[j]*(b*b/2.0*(thisQ[j]-deltaQ*a/(b-a))+b*b*b*deltaQ/3.0/(b-a)
 //            -                 a*a/2.0*(thisQ[j]-deltaQ*a/(b-a))-a*a*a*deltaQ/3.0/(b-a));

//            g = log(F0half[j-1]/F0half[j+1])/(Egrid.Ecc[j+1]-Egrid.Ecc[j-1]);           
//            c = thisQ[j]-a*deltaQ/(b-a);
//            d = deltaQ/(b-a); 
//            IntExp(exp1,a,b,g,1);
//            IntExp(exp2,a,b,g,2);          
//            PU = Ng*gamma*F0[j]*exp(Egrid.Ecc[j]*g)*(c*exp1+d*exp2); // method 3
          
            // Update RHS dint(f)/dE/dt due to source term
            //  
            ExcS[j] -= PU; //Egrid.dE;
            ExcS[0] += PU; //Egrid.dE;
            deltaE = thisU - Egrid.Ece[j];
            //cout << "this J = " << j << endl;
         }
         if(Egrid.Ece[j] >= thisU) {
            //cout << "this j = " << j << endl;
            Esub = Egrid.Ece[j]+deltaE;
  //          deltaQ = thisQ[j+1]-thisQ[j];
            
            // calculate source term for cell j
            //
            a = Egrid.Ece[j];
            b = Esub;
            PL = Ng*gamma*F0half[j]*(thisQ[j]+thisQ[j+1])/2.0*(b*b-a*a)/2.0; // method 1
           // PL = Ng*gamma*F0half[j]*(b*b/2.0*(thisQ[j]-deltaQ*a/Egrid.dE)+b*b*b*deltaQ/3.0/Egrid.dE
           //    -                 a*a/2.0*(thisQ[j]-deltaQ*a/Egrid.dE)-a*a*a*deltaQ/3.0/Egrid.dE);
            a = Esub;
            b = Egrid.Ece[j+1];
            PU = Ng*gamma*F0half[j]*(thisQ[j]+thisQ[j+1])/2.0*(b*b-a*a)/2.0; // method 1
           // PU = Ng*gamma*F0half[j]*(b*b/2.0*(thisQ[j]-deltaQ*a/Egrid.dE)+b*b*b*deltaQ/3.0/Egrid.dE
           //    -                 a*a/2.0*(thisQ[j]-deltaQ*a/Egrid.dE)-a*a*a*deltaQ/3.0/Egrid.dE);

            // begin method 3 (See Hagelaar Eqs 46-51 (note sign error in 51))
            //
//             g = log(F0half[j-1]/F0half[j+1])/(Egrid.Ecc[j+1]-Egrid.Ecc[j-1]);
//             if(j==nE-1) g = log(F0half[j-1]/F0half[j])/(Egrid.Ecc[j+1]-Egrid.Ecc[j]);
//             a = Egrid.Ece[j];
//             b = Esub;
//             c = thisQ[j]-a*deltaQ/Egrid.dE;
//             d = deltaQ/Egrid.dE; 
//             IntExp(exp1,a,b,g,1);
//             IntExp(exp2,a,b,g,2);          
//             PL = Ng*gamma*F0[j]*exp(Egrid.Ecc[j]*g)*(c*exp1+d*exp2);
//             
//             a = Esub;
//             b = Egrid.Ece[j+1];
//             c = thisQ[j]-a*deltaQ/Egrid.dE;
//             d = deltaQ/Egrid.dE; 
//             IntExp(exp1,a,b,g,1);
//             IntExp(exp2,a,b,g,2); 
//             PU = Ng*gamma*F0[j]*exp(Egrid.Ecc[j]*g)*(c*exp1+d*exp2);
            //
            // end method 3

            // Update RHS dint(f)/dE/dt due to excitation source term
            //  
            ExcS[j] -= (PL+PU); 
            ExcS[thism-1] += PL; 
            ExcS[thism] += PU;
            thism = thism+1;
         } 
      }
   }
   for (auto j=0; j<nE; j++) {
      ExcS[j] /= Egrid.dE*sqrt(Egrid.Ecc[j]);
   }

}



void EEDF::advanceF0(const energyGrid& Egrid, const double& dt)
{
   const int nE = Egrid.nE;
   /*
   // Explicit forward advance
   //
   for (auto n=0; n<nE; n++) {
      F0[n] = F0[n] - dt*(Flux[n+1]-Flux[n])/Egrid.dE/sqrt(Egrid.Ecc[n])
                    + dt*ExcS;
   }
   */
   
   // Implicit Crank-Nicolson advance
   //
   
   // Compute the tridiagonal matrix components
   //
   vector<double> a(nE,0.0), b(nE,0.0), c(nE,0.0), d(nE,0.0);
   setTriDiagCoeffs(a,b,c,d,Egrid,dt);
   
   // add inelastic source/sink term to d vector
   //
   for (auto n=0; n<nE; n++) {
      d[n] += 4.0*dt*ExcS[n];
   }
   
   // Row reduce
   //
   for (auto n=1; n<nE; n++) {
      b[n] = b[n] - a[n]*c[n-1]/b[n-1];
      d[n] = d[n] - a[n]*d[n-1]/b[n-1];
   }
   
   // Back solve
   //
   F0[nE-1] = d[nE-1]/b[nE-1];
   for (auto n=nE-2; n>-1; n--) {
      F0[n] = (d[n] - F0[n+1]*c[n])/b[n]; 
      if(F0[n]<0) {
         cout << "F0("<< n <<")= " << F0[n] << endl;
         F0[n] = 1.0e-15;
      }
   }
   
   // update moments and F0half
   //
   Te0 = 0;
   zeroMom = 0;
   for (auto n=0; n<nE; n++) {
      F0half[n] = (F0[n]+F0old[n])/2.0;
      Te0 = Te0 + 2.0/3.0*pow(Egrid.Ecc[n],1.5)*F0[n]*Egrid.dE;
      zeroMom = zeroMom + sqrt(Egrid.Ecc[n])*F0[n]*Egrid.dE;
   }
   
}

void EEDF::IntExp(double& soln, const double& a, const double& b, const double& g, const int& pow)
{ 
   if(pow==1) {
      soln = -( exp(-b*g)*(b*g+1.0) - exp(-a*g)*(a*g+1.0) )/(g*g);
   }
   if(pow==2) {
      soln = -( exp(-b*g)*(b*b*g*g + 2.0*(b*g+1.0) ) 
           -    exp(-a*g)*(a*a*g*g + 2.0*(a*g+1.0) ) )/(g*g*g);
   }
   if(pow!=1 && pow!=2) {
      cout << "Error: call to IntExp only accepts pow=1 and 2" << endl;
      exit(EXIT_FAILURE);
   }

}
void EEDF::setTriDiagCoeffs(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d,
                            const energyGrid& Egrid, const double& dt)
{  
   // create vectors of Courant and alphas
   //
   const int nE = Egrid.nE;
   vector<double> couraM(nE,0.0), alphaM(nE,0.0); // W and D at n-1/2
   vector<double> couraP(nE,0.0), alphaP(nE,0.0); // W and D at n+1/2
      
   for(auto n=0; n<nE; n++) {
      couraM[n] = W[n]*dt/Egrid.dE/sqrt(Egrid.Ecc[n]);
      alphaM[n] = D[n]*dt/(Egrid.dE*Egrid.dE*sqrt(Egrid.Ecc[n]));
      couraP[n] = W[n+1]*dt/Egrid.dE/sqrt(Egrid.Ecc[n]);
      alphaP[n] = D[n+1]*dt/(Egrid.dE*Egrid.dE*sqrt(Egrid.Ecc[n]));  
   }
   
   // set tridiag coefficients 
   //
   for(auto n=0; n<nE; n++) {
      a[n] = - (couraM[n] + 2.0*alphaM[n]);
      b[n] = 4.0 + couraP[n] - couraM[n] + 2.0*(alphaP[n] + alphaM[n]);
      c[n] = couraP[n] - 2.0*alphaP[n];     
   }
   //const double EmaxBC = 0; // zero value at Emax+1/2
   const double EmaxBC = (2.0*D[nE]+Egrid.dE*W[nE])/(2.0*D[nE]-Egrid.dE*W[nE]); // 0 Flux at Emax
   b[nE-1] += c[nE-1]*EmaxBC;
      
   // set d coefficient A*F0^{n+1} = d(F0^{n})
   //
   d[0] = (4.0 - couraP[0] - 2.0*alphaP[0])*F0old[0] 
        - (couraP[0] - 2.0*alphaP[0])*F0old[1];
   for(auto n=1; n<nE-1; n++) {
      d[n] = (couraM[n] + 2*alphaM[n])*F0old[n-1]
           + (4.0 + couraM[n] - couraP[n] - 2.0*(alphaM[n] + alphaP[n]))*F0old[n]
           - (couraP[n] - 2.0*alphaP[n])*F0old[n+1];
   }
   d[nE-1] = (couraM[nE-1]+2.0*alphaM[nE-1])*F0old[nE-2]
           + (4.0+couraM[nE-1]-couraP[nE-1]-2.0*(alphaM[nE-1]+alphaP[nE-1]))*F0old[nE-1];
   d[nE-1] -= (couraP[nE-1]-2.0*alphaP[nE-1])*EmaxBC*F0old[nE-1];
}

#endif
