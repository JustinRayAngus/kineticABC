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
#include "electricField.h"
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
  
 
   // set electric field E[V/m]
   //
   electricField ElcField;
   ElcField.initialize(inputRoot, dataFile);
   const double EVpm = ElcField.EVpm; 

   // initialize eedf
   //  
   EEDF eedf;
   eedf.initialize(gas, Egrid, inputRoot, dataFile);   

   // compute flux at cell-edges using initial F=F0old=F0half
   //
   eedf.computeFlux(gas, Egrid, EVpm);
   dataFile.add(eedf.Flux, "Flux", 1);
   tDom.setdtSim(eedf, Egrid); // set initial time step


   // march forward in time
   //
   double dtSim = tDom.dtSim;
   cout << "Initial simulation time step: " << dtSim << endl << endl;
   double thist = 0;
   int thistOutInt = 1;
   vector<double> F0m(Egrid.nE,0.0), errorVec(Egrid.nE,0.0);
   double error;

   while(thist<tDom.tmax) {
      thist = thist + dtSim; // new time at end of this time step
      for (auto i=0; i<100; i++) {     // predict-correct iterations
         eedf.advanceF0(Egrid, dtSim); // F0 and F0half updated here
         if (i==0) { 
            F0m = eedf.F0;
         }
         else {
            for (auto j=0; j<Egrid.nE; j++) {
               errorVec[j] = abs(1.0-F0m[j]/eedf.F0[j]);
            }
            error = *max_element(begin(errorVec), end(errorVec));
            if (error <= 1e-5) {
               if (i>=10) {
                  cout << "iteration = " << i << endl;
                  cout << "error = " << error << endl;
               }
               break;
            }
            F0m = eedf.F0;
         }
         eedf.computeIznS(gas, Egrid); // uses F0half
         eedf.computeFlux(gas, Egrid, EVpm);
         eedf.computeExcS(gas, Egrid); // uses F0half
      }

      // check if thist is an output time
      //
      if(thist >= tDom.tOutVec[thistOutInt]) {
         tDom.updatetOut(thist);
         dataFile.writeAll(); // append extendable outputs
         cout << "Output variables dumped at t = " << thist << " s" << endl;
         thistOutInt = thistOutInt+1;
         //cout << "Simulation time step = " << dtSim << endl;
      }

      // reset F0old, F0half, and update simulation time step
      //
      eedf.F0old = eedf.F0;
      eedf.F0half = eedf.F0;
      eedf.computeIznS(gas, Egrid);
      eedf.computeFlux(gas, Egrid, EVpm);
      eedf.computeExcS(gas, Egrid);
      tDom.setdtSim(eedf, Egrid);
      dtSim = tDom.dtSim;
   }

   cout << "\nFinal simulation time step = " << dtSim << endl;
   cout << "\nEnding simulation" << endl;
   cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" << endl;
   return 1;


   //testing::InitGoogleMock(&argc, argv);
   //return RUN_ALL_TESTS();
}


