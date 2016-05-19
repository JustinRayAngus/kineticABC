#include "gmock/gmock.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <typeinfo>
#include "json/json.h"
#include "H5Cpp.h"

#include "energyGrid.h"
#include "timeDomain.h"
#include "Gas.h"
#include "HDF5dataFile.h"
#include "EEDF.h"

using namespace std;

#ifndef H5_NO_NAMESPACE
   using namespace H5;
#endif


int main(int argc, char** argv) {   
   cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
   cout << "" << endl;
   cout << "Initiating simulation" << endl;

   
   // Parse the specified json input file
   //
   const string inputFile = "../input.json";
   Json::Value inputRoot; // will contain root value after parsing
   Json::Reader reader; 
   ifstream ifile(inputFile);
   bool isJsonOK = (ifile !=NULL && reader.parse(ifile, inputRoot));
   if(isJsonOK) {
      printf("\nInput file %s %s", inputFile.c_str(), "parsed successfully \n");
   }
   else {
      printf("ERROR: json input file not found or cannot be parsed due to errors \n");
      exit (EXIT_FAILURE);
   }


   // set output file in HD5FdataFile class
   //
   HDF5dataFile dataFile;
   const string outputFile = "output.h5";
   dataFile.setOutputFile(outputFile);

   
   // initialize energy grid, time domain, and EEDF
   //
   energyGrid Egrid;
   Egrid.initialize(inputRoot, dataFile);
   //
   timeDomain tDom;
   tDom.initialize(inputRoot, dataFile);
   //
   Gas gas;
   gas.initialize(Egrid, inputRoot, dataFile);   
   //
   EEDF eedf;
   eedf.initialize(Egrid, inputRoot, dataFile);   

   
   // determine E[V/m] for const Qelm
   //
   const double EVpm = eedf.Te0*gas.Ng*gas.Qelm[1]*sqrt(3.0*gas.mM); 
   cout << "E = " << EVpm << " [V/m] " << endl << endl;
     

   // compute flux at cell-edges using initial F0
   //
   eedf.computeFlux(gas, Egrid, EVpm);


   // write outputs at specified times
   //
   const double dtSim = eedf.dtStable*1e3;
   cout << "Stable time step: " << eedf.dtStable << endl;
   cout << "Simulation time step: " << dtSim << endl << endl;        
   double thist = 0;
   int thistOutInt = 1;

   while(thist<tDom.tmax) {
      thist = thist + dtSim;
      eedf.advanceF0(Egrid, dtSim);
      eedf.computeFlux(gas, Egrid, EVpm);
      if(thist >= tDom.tOutVec[thistOutInt]) {
         tDom.updatetOut(thist);
         dataFile.writeAll(); // append extendable outputs
         cout << "Output variables dumped at time " << thist << endl;
         thistOutInt = thistOutInt+1;
      }
   }
   //dataFile.writeAll(); // appends extendable outputs
   //dataFile.writeAll(); // appends extendable outputs again


   cout << "\nEnding simulation" << endl;
   cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" << endl;
   return 1;


   //testing::InitGoogleMock(&argc, argv);
   //return RUN_ALL_TESTS();
}


