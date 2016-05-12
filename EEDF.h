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

using namespace std;

class EEDF
{
public:
   double Te0;        // initial Te0;
   string type0;      // initial type of EEDF
   vector<double> f0; // initial EEDF at cell-center 
   void initialize(const energyGrid&, const string&);
};

void EEDF::initialize(const energyGrid& Egrid, const string& inputFile)
{
   Json::Value root; // will contain root value after parsing
   Json::Reader reader;
   const Json::Value defValue; // used for default reference
 
   ifstream ifile(inputFile);
   bool isJsonOK = (ifile !=NULL && reader.parse(ifile, root) );
   if(isJsonOK) {
      const Json::Value EEDF = root.get("EEDF",defValue);
      if(EEDF.isObject()) {
         printf("Reading EEDF information from file: %s %s",
                inputFile.c_str(), "\n");
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
         cout << "initial EEDF is " << type0 << " with Te = " << Te0 << endl << endl;
      }
      else {
         cout << "value for key \"Egrid\" is not object type !" << endl;
      }
   }
   else {
      printf("ERROR: json is not OK ... input file not found? \n");
      exit (EXIT_FAILURE);
   }
   f0.assign(Egrid.nE,0.0);
   //f0ce.assign(Egrid.nE+1,0.0);
   for (auto n=0; n<Egrid.nE; n++) {
      f0[n] = 2.0/sqrt(Te0*3.14159)/Te0/Te0*exp(-Egrid.Ecc[n]/Te0);
      //f0ce[n] = n*dE;
   }
   //f0ce[Egrid.nE] = Emax;
   // cout << Ecc.size() << endl; // always gives 0????
}


#endif
