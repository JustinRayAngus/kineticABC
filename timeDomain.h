/***
 * 
 * time domain class
 *
***/

#ifndef timeDomain_h
#define timeDomain_h

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "json/json.h"
#include "HDF5dataFile.h"

using namespace std;

class timeDomain
{
 
public:
  double dtOut, tOutSteps;   // Output intervals and number of steps
  double Cm, tmax;           // Courant multiplier and max time
  vector<double> tOutVec;        // vector of output times
  double tOut;           // current output time 
  void initialize(const Json::Value&, HDF5dataFile&);
  
  void updatetOut(double thistOut);

};

void timeDomain::initialize(const Json::Value& root, HDF5dataFile& dataFile)
{
   const Json::Value defValue; // used for default reference
   const Json::Value Time = root.get("Time",defValue);
   if(Time.isObject()) {
      printf("Initializing time domain ...\n");
      Json::Value dtOutVal      = Time.get("dtOut",defValue);
      Json::Value tOutStepsVal  = Time.get("tOutSteps",defValue);
      Json::Value CmVal         = Time.get("Cm",defValue);
      if(dtOutVal == defValue || CmVal == defValue || tOutStepsVal == defValue) {
         printf("ERROR: default 'Time' variables not declared in input file\n");
         exit (EXIT_FAILURE);
      }
      dtOut = dtOutVal.asDouble();
      tOutSteps = tOutStepsVal.asDouble();
      Cm = CmVal.asDouble();
      tmax = dtOut*tOutSteps;
      if(Cm >= 1 ) {
         printf("WARNING: Courant number is >= 1 !!! \n");
         //exit (EXIT_FAILURE);
      }
      //
      cout << "tmax = " << dtOut*tOutSteps << endl;
      cout << "tOut intervals = " << dtOut << endl;
      cout << "Courant multiplier = " << Cm << endl;
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



#endif
