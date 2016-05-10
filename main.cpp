#include "gmock/gmock.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <typeinfo>
#include "json/json.h"
#include "energyGrid.h"
#include "HDF5FileTools.h"
#include "H5Cpp.h"

using namespace std;

#ifndef H5_NO_NAMESPACE
   using namespace H5;
#endif

const int NX = 5;
const int NY = 6;
const int RANK = 2;

int main(int argc, char** argv) {
   
   cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
   cout << "" << endl;
   cout << "Initiating simulation" << endl;
   
   // Create energy grid from a specified json input file
   //
   string inputFile = "input_file.json";
   energyGrid Egrid;
   Egrid.setEnergyGrid(inputFile);
   cout << "nE = " << Egrid.nE << endl;
   cout << "Emax = " << Egrid.Emax << endl;
   cout << "dE = " << Egrid.dE << endl << endl; 
   //for (auto n = 0; n < Egrid.nE; n++) {
   //   cout << "Ecc[" << n << "] = " << Egrid.Ecc[n] << endl; 
   //   cout << "Ece[" << n << "] = " << Egrid.Ece[n] << endl; 
   //}

   
   // load more information and manipulate f0 in time
   //
   vector<vector<double>> vec2d;
   int nrow = 5;
   int ncol = 4;
   vec2d.resize(nrow,vector<double>(ncol,0.0));
   //cout << vec2d.size() << endl;
   //cout << vec2d[0].size() << endl;


   // write energy grid to a specifid h5 output file
   // 
   HDF5FileTools hdf5FileTools;
   const string outputFile = "output.h5";
   H5File file( outputFile.c_str(), H5F_ACC_TRUNC ); // overwrites old
   //
   const string varNameEcc = "Ecc";
   hdf5FileTools.writeOutput(outputFile, varNameEcc, Egrid.Ecc);
   //
   const string varNameEce = "Ece";
   hdf5FileTools.writeOutput(outputFile, varNameEce, Egrid.Ece);
   //
   //const string varNamevec2d = "vec2d";
   //hdf5FileTools.writeOutput(outputFile, varNamevec2d, vec2d);
   

   cout << "\nEnding simulation" << endl;
   cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
   cout << "" << endl;
   return 1;

   //testing::InitGoogleMock(&argc, argv);
   //return RUN_ALL_TESTS();
}


