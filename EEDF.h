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
#include "json/json.h"
#include "energyGrid.h"
#include "HDF5dataFile.h"

using namespace std;

class EEDF
{
public:
   double Te0;        // initial Te0;
   string type0;      // initial type of EEDF
   vector<double> f0; // initial EEDF at cell-center 
   void initialize(const energyGrid&, const Json::Value&, HDF5dataFile&);
};

void EEDF::initialize(const energyGrid& Egrid, const Json::Value& root, HDF5dataFile& dataFile)
{
   f0.assign(Egrid.nE,0.0);
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
      if(type0=="gaussian" || type0=="Gaussian") {
         cout << "Initial EEDF is Gaussian with Te = " << Te0 << endl;
         for (auto n=0; n<Egrid.nE; n++) {
            f0[n] = 2.0/sqrt(Te0*3.14159)/Te0*exp(-Egrid.Ecc[n]/Te0);
            //f0ce[n] = n*dE;
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
   
   dataFile.add(f0, "f0", 1);
   dataFile.add(Te0, "Te0", 1); 
   cout << endl;  
}


#endif
