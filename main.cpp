#include "gmock/gmock.h"
#include <iostream>
#include <string>
#include <fstream>
#include "json/json.h"
#include "energyGrid.h"

using namespace std;

int main(int argc, char** argv) {

   // Read energy grid from specified input file
   //
   string inputFile = "input_file.json";
   energyGrid Egrid;
   Egrid.setEnergyGrid(inputFile);
   cout << "nE = " << Egrid.nE << endl;
   cout << "Emax = " << Egrid.Emax << endl;
   cout << "dE = " << Egrid.dE << endl; 
   //for (auto n = 0; n < Egrid.nE; n++) {
   //   cout << "Ecc[" << n << "] = " << Egrid.Ecc[n] << endl; 
   //   cout << "Ece[" << n << "] = " << Egrid.Ece[n] << endl; 
   //}
   
   // now load other information and manipulate f0 in time
   //

   // now write stuff to hdf5 files
   //
   cout << "\n next step is to write stuff to hdf5 file" << endl;   

   return 1;

   //testing::InitGoogleMock(&argc, argv);
   //return RUN_ALL_TESTS();
}

