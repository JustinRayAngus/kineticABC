/***
 * 
 * energy grid class
 *
***/

#ifndef energyGrid_h
#define energyGrid_h

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "json/json.h"

using namespace std;

class energyGrid
{
public:
   int nE;
   double Emax, dE;
   vector<double> Ecc, Ece; // Energy at cell-center and at cell-edge 
   //void setEnergyGrid(const string&);
   void initialize(const Json::Value&);
};

void energyGrid::initialize(const Json::Value& root)
{
   const Json::Value defValue; // used for default reference
   const Json::Value Egrid = root.get("Egrid",defValue);
   if(Egrid.isObject()) {
      printf("\nInitializing energy grid ...\n");
      Json::Value EmaxVal = Egrid.get("Emax",defValue);
      Json::Value nEVal = Egrid.get("nE",defValue);
      if(EmaxVal == defValue || nEVal == defValue) {
         printf("ERROR: Emax or nE not declared in input file\n");
         exit (EXIT_FAILURE);
      } 
      nE = nEVal.asInt();
      if(nE != nEVal.asDouble() || nE < 1) {
         printf("ERROR: nE is not set as a positive integer in input file\n");
         exit (EXIT_FAILURE);
      }
      Emax = EmaxVal.asDouble();
      dE = Emax/nE;
      //
      cout << "nE = " << nE << endl;
      cout << "Emax = " << Emax << endl;
      cout << "dE = " << dE << endl << endl;
      //for (auto n = 0; n < Egrid.nE; n++) {
      //   cout << "Ecc[" << n << "] = " << Egrid.Ecc[n] << endl;
      //   cout << "Ece[" << n << "] = " << Egrid.Ece[n] << endl;
      //}
   }
   else {
      cout << "value for key \"Egrid\" is not object type !" << endl;
   }
  
   //Ecc.reserve(nE);
   //Ece.reserve(nE+1);
   Ecc.assign(nE,0.0);
   Ece.assign(nE+1,0.0);
   for (auto n=0; n<nE; n++) {
      Ecc[n] = 0.5*dE+n*dE;
      Ece[n] = n*dE;
   }
   Ece[nE] = Emax;
   // cout << Ecc.size() << endl; // always gives 0????
}


#endif
