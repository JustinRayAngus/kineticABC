/***
 * 
 * time domain class
 *
***/

#ifndef timeDomain_h
#define timeDomain_h

#include "EEDF.h"
#include "energyGrid.h"

using namespace std;

class timeDomain
{
 
public:
  double dtOut, tOutSteps;   // Output intervals and number of steps
  double tmax;           // Max time
  double dtSim;              // Simulation time-step
  double dtFrac;             // dtSim = dtmax/dtFrac 
  vector<double> tOutVec;        // vector of output times
  double tOut;           // current output time 
  
  void initialize(const Json::Value&, HDF5dataFile&);
  void updatetOut(double thistOut);
  void setdtSim(const EEDF&, const energyGrid&); // set time-step for implicit solver

};

void timeDomain::initialize(const Json::Value& root, HDF5dataFile& dataFile)
{
   const Json::Value defValue; // used for default reference
   const Json::Value Time = root.get("Time",defValue);
   if(Time.isObject()) {
      printf("Initializing time domain ...\n");
      Json::Value dtOutVal      = Time.get("dtOut",defValue);
      Json::Value tOutStepsVal  = Time.get("tOutSteps",defValue);
      Json::Value dtFracVal     = Time.get("dtFrac",defValue);
      if(dtOutVal == defValue || tOutStepsVal == defValue || dtFracVal == defValue) {
         printf("ERROR: default 'Time' variables not declared in input file\n");
         exit (EXIT_FAILURE);
      }
      dtOut = dtOutVal.asDouble();
      tOutSteps = tOutStepsVal.asDouble();
      dtFrac = dtFracVal.asDouble();
      tmax = dtOut*tOutSteps;
      //
      cout << "tmax = " << dtOut*tOutSteps << endl;
      cout << "tOut intervals = " << dtOut << endl;
      cout << "dtFrac = " << dtFrac << endl;
   }
   else {
      cout << "value for key \"Time\" is not object type !" << endl;
   }
  
   tOutVec.assign(tOutSteps+1,0.0);
   for(auto n=0; n<tOutSteps+1; n++) {
      tOutVec[n] = n*dtOut;  // vector of output times for references
      //cout << tOutVec[n] << endl;
   }
   
   tOut = 0.0; // first output is always at t=0
   dataFile.add(tOut, "tout", 1); // actual output time
   cout << endl;
}

void timeDomain::updatetOut(const double thistOut)
{
   tOut = thistOut;
}

void timeDomain::setdtSim(const EEDF& eedf, const energyGrid& Egrid)
{
   const vector<double> ExcS = eedf.ExcS;
   const vector<double> IznS = eedf.IznS;
   double nunet = eedf.nunet;
   const vector<double> Flux = eedf.Flux;
   const vector<double> Ecc = Egrid.Ecc;
   
   // F0^{n+1} = F0^n + dt*(S - divFlux)
   //
   const int nE = Egrid.nE;
   double thisRHS;
   double dtmax = 1.0e6;
   for(auto i=0; i<nE; i++) {
      thisRHS = ExcS[i] + IznS[i] - nunet*eedf.F0half[i]
             - (Flux[i+1]-Flux[i])/(sqrt(Ecc[i])*Egrid.dE);
      if(thisRHS!=0) {
         dtmax = min(dtmax, abs(eedf.F0old[i]/thisRHS));
      }
   }
   dtSim = min(dtOut/dtFrac,dtmax/dtFrac);
   //cout << "dtSim = " << dtSim << endl;
}

#endif
