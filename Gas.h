/***
 * 
 * background Gas class
 *
***/

#ifndef Gas_h
#define Gas_h

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include "json/json.h"
#include "energyGrid.h"
#include "HDF5dataFile.h"

using namespace std;

class Gas
{
public:
   double Tg, Pg, Ng;   // gas temp [K], press [Torr], and dens [1/m^3]
   double mM;           // electron mass / gas species mass
   string gasName;      // gas species name
   vector<double> Qelm; // elastic-momentum-transfer cross section [m^2] 
                        // at cell edges
                        
   void initialize(const energyGrid&, const Json::Value&, HDF5dataFile&);
};

void Gas::initialize(const energyGrid& Egrid, const Json::Value& root, HDF5dataFile& dataFile)
{
   Qelm.assign(Egrid.nE+1,0.0);
   const Json::Value defValue; // used for default reference
   const Json::Value Gas = root.get("Gas",defValue);
   if(Gas.isObject()) {
      printf("Initializing background gas ...\n");
      Json::Value TgVal = Gas.get("Tg",defValue);
      Json::Value PgVal = Gas.get("Pg",defValue);
      Json::Value nameVal = Gas.get("name",defValue);
      if(TgVal == defValue || PgVal == defValue || nameVal == defValue) {
         printf("ERROR: Tg, Pg, or name of gas not specified correctly in input file\n");
         exit (EXIT_FAILURE);
      } 
      Tg = TgVal.asDouble(); // [K]
      Pg = PgVal.asDouble();
      Ng = Pg/760*273/Tg*2.6868e25; // gas density [1/m^3]
      if(Tg < 0.0) {
         printf("ERROR: Tg is not a positive value in input file\n");
         exit (EXIT_FAILURE);
      }
      gasName = nameVal.asString();
      if(gasName=="nitrogen" || gasName=="Nitrogen") {
      
         const double atomicM = 14.0067; // atomic mass of nitrogen atom
         const double meconst = 9.1094e-28;   // electron mass [g]
         const double amuconst = 1.6605e-24;  // atomic mass unit [g]
         mM = meconst/(amuconst*atomicM*2.0); // me/MN2
         
         cout << "Background gas is nitrogen ... " << endl;
         cout << "Tg = " << Tg << " K" << endl;
         cout << "Pg = " << Pg << " Torr" << endl;
         cout << "Ng = " << Ng << " 1/m^3" << endl;
         cout << "m/M = " << mM << endl;
         const int nQelm = Qelm.size();
         for (int n=0; n<nQelm; n++) {
            Qelm[n] = 1e-20;
         }
      }
      else {
         cout << "Gas species " << gasName << " is not a valid option" << endl;
         exit (EXIT_FAILURE);
      }
   }
   else {
      cout << "value for key \"Gas\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   
   dataFile.add(Tg, "Tg", 0);
   dataFile.add(Pg, "Pg", 0);
   dataFile.add(Ng, "Ng", 0);
   dataFile.add(mM, "mM", 0);
   dataFile.add(Qelm, "Qelm", 0);
    
   cout << endl;  
}


#endif
