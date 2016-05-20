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
   vector<double> F0;     // EEDF
   vector<double> W, D, Flux;  // Energy space adv, diff, and flux at cell-edge
   
   double dtStable;  // stable time-step
   
   void initialize(const energyGrid&, const Json::Value&, HDF5dataFile&);
   void computeFlux(const Gas&, const energyGrid&, const double&);
   void advanceF0(const energyGrid&, const double&);
};

void EEDF::initialize(const energyGrid& Egrid, const Json::Value& root, HDF5dataFile& dataFile)
{
   F0.assign(Egrid.nE,0.0);
   Flux.assign(Egrid.nE+1,0.0);
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
   
   dataFile.add(F0, "F0", 1);  
   dataFile.add(Flux, "Flux", 1); //Flux = WF - D*dF/dE
   dataFile.add(W, "W", 1);  // advection coefficient
   dataFile.add(D, "D", 1);  // diffusion coefficient
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
   const double gamma = sqrt(2*econst/meconst);
   const int nE = Egrid.nE; // number of cell-center points
   
   // compute advection and diffusion coefficients
   //
   vector<double> PecNum, dtMaxW, dtMaxD, dtMaxWD;
   PecNum.assign(nE+1,0.0);    // W*dE/D, Peclet number
   dtMaxW.assign(nE+1,1E20);    // Courant: C=W*dt/dE < 1
   dtMaxD.assign(nE+1,1E20);    // Diff:    D=U*dt/dE^2 < 1/2
   dtMaxWD.assign(nE+1,1E20);  // C^2/D < 2
   
   for (auto n=1; n<nE+1; n++) { // dont need values at E=0
      W[n] = -gamma*Ng*2.0*mM*pow(Egrid.Ece[n],2.0)*gas.Qelm[n];
      D[n] = gamma*EVpm*EVpm/3.0/Ng/gas.Qelm[n]*Egrid.Ece[n];
      PecNum[n] = W[n]*Egrid.dE/D[n];
      dtMaxW[n] = -Egrid.dE/W[n];
      dtMaxD[n] = 0.5*Egrid.dE*Egrid.dE/D[n];
      dtMaxWD[n] = 2.0*abs(D[n])/W[n]/W[n];
   }
   double dtW = *min_element(begin(dtMaxW), end(dtMaxW));
   double dtD = *min_element(begin(dtMaxD), end(dtMaxD));
   double dtWD = *min_element(begin(dtMaxWD), end(dtMaxWD));
   //cout << "Stable time step for advection: " << dtW << endl;
   //cout << "Stable time step for diffusion: " << dtD << endl;
   //cout << "Stable time step for advection-diffusion: " << dtWD << endl;
   
   dtStable = min(min(dtW,dtD),dtWD);
   //cout << "dtStable = " << dtStable*1e9 << " ns" << endl;
   //cout << endl;
   
   //cout << PecNum[0] << endl;
   
   // compute gradient of F0
   //
   vector<double> DF0DE(nE+1,0.0);  // B.C. : dF0/dE = 0 at E=0
   vector<double> F0_ce(nE+1,0.0);  // F0 at cell-edge
   for (auto n=1; n<nE; n++) {
      DF0DE[n] = (F0[n]-F0[n-1])/Egrid.dE;
      F0_ce[n] = (F0[n]+F0[n-1])/2.0;
   }
   DF0DE[nE] = (0.0 - F0[nE-1])/Egrid.dE; // B.C.: F0 = 0 at Emax+1/2*dE
   F0_ce[nE] = F0[nE-1]/2.0;
      
   // compute Flux = D*dF0/dE - W*F0
   //
//    for (auto n=1; n<nE; n++) {
//       Flux[n] = W[n]*F0[n-1]/(1.0-exp(-PecNum[n])) + W[n]*F0[n]/(1.0-exp(PecNum[n])); 
//    }
//    Flux[nE] = W[nE]*F0[nE-1]/(1.0-exp(-PecNum[nE]));

   for (auto n=1; n<nE; n++) { 
      Flux[n] = W[n]*F0_ce[n] - D[n]*DF0DE[n];
   }
   Flux[nE] = W[nE]*F0_ce[nE] - D[nE]*DF0DE[nE];
   
   
}

void EEDF::advanceF0(const energyGrid& Egrid, const double& dt)
{
   Te0 = 0;
   zeroMom = 0;
   const int nmax= F0.size();
   for (auto n=0; n<nmax; n++) {
      F0[n] = F0[n] + dt*(Flux[n]-Flux[n+1])/Egrid.dE/sqrt(Egrid.Ecc[n]);
      Te0 = Te0 + 2.0/3.0*pow(Egrid.Ecc[n],1.5)*F0[n]*Egrid.dE;
      zeroMom = zeroMom + sqrt(Egrid.Ecc[n])*F0[n]*Egrid.dE;
   }
   
}


#endif
