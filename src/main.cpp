#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <typeinfo>
#include <algorithm>

#include "json/json.h"
#include "HDF5dataFile.h"
#include "energyGrid.h"
#include "timeDomain.h"
#include "Gas.h"
#include "electricField.h"
#include "EEDF.h"

using namespace std;

#ifndef H5_NO_NAMESPACE
   using namespace H5;
#endif


int main(int argc, char** argv) {   
   cout << endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
   cout << "" << endl;
   cout << "Initiating simulation" << endl;


   // Parse the specified json input file
   //
   const string inputFile = "./input.json";
   Json::Value inputRoot; // will contain root value after parsing
   Json::Reader reader; 
   ifstream ifile(inputFile);
   bool isJsonOK = (ifile !=NULL && reader.parse(ifile, inputRoot));
   if(isJsonOK) {
      cout << "Input " << inputFile << " parsed successfully" << endl;   
   }
   else {
      cout << "ERROR: json input file not found or cannot be parsed due to errors" << endl;
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
   Egrid.initialize(inputRoot);
   dataFile.add(Egrid.Ecc, "Ecc", 0);
   dataFile.add(Egrid.Ece, "Ece", 0); 
  
   //
   timeDomain tDom;
   tDom.initialize(inputRoot);
   dataFile.add(tDom.tOut, "tout", 1); // actual output time   

   //
   Gas gas;
   gas.initialize(Egrid, inputRoot);   
   dataFile.add(gas.Tg,   "Tg", 0); 
   dataFile.add(gas.Pg,   "Pg", 0); 
   dataFile.add(gas.Ng,   "Ng", 0); 
   dataFile.add(gas.Qelm, "Qelm", 0); 
   dataFile.add(gas.mM,   "mM", 0); 
   dataFile.add(gas.Qexc, "Qexc", 0); 
   dataFile.add(gas.Uexc, "Uexc", 0); 
   dataFile.add(gas.Qizn, "Qizn", 0); 
   dataFile.add(gas.Uizn, "Uizn", 0); 
   dataFile.add(gas.Qmom, "Qmom", 0); 

   //
   electricField ElcField;
   ElcField.initialize(inputRoot);
   dataFile.add(ElcField.EVpm, "E", 0);

   // initialize eedf
   //  
   EEDF eedf;
   eedf.initialize(gas, Egrid, inputRoot);   
   dataFile.add(eedf.F0, "F0", 1);
   dataFile.add(eedf.ExcS, "ExcS", 1); // excitation source term
   dataFile.add(eedf.IznS, "IznS", 1); // ionization source term
   dataFile.add(eedf.W, "W", 1);       // advection coefficient
   dataFile.add(eedf.D, "D", 1);       // diffusion coefficient
   dataFile.add(eedf.zeroMom, "zeroMom", 1); // should remain at one always !!! 
   dataFile.add(eedf.Te, "Te", 1);
   dataFile.add(eedf.nunet, "nunet", 1); // dne/dt/ne (electron production rate)


   // compute flux at cell-edges using initial F=F0old=F0half
   //
   eedf.computeFlux(gas, Egrid, ElcField.EVpm);
   dataFile.add(eedf.Flux, "Flux", 1);
   cout << endl;
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
            //if (i==6) {
               if (i>=10) {
                  cout << "iteration = " << i << endl;
                  cout << "error = " << error << endl;
               }
               break;
            }
            F0m = eedf.F0;
         }
         eedf.computeIznS(gas, Egrid); // uses F0half
         eedf.computeExcS(gas, Egrid); // uses F0half
         eedf.computeFlux(gas, Egrid, ElcField.EVpm);
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
      eedf.computeExcS(gas, Egrid); // uses F0half
      eedf.computeFlux(gas, Egrid, ElcField.EVpm);
      tDom.setdtSim(eedf, Egrid);
      dtSim = tDom.dtSim;
   }

   cout << endl << "Final simulation time step = " << dtSim << endl;
   cout << endl << "Ending simulation" << endl;
   cout << endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl;
   return 1;
}


