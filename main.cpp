#include "gmock/gmock.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <typeinfo>
#include "json/json.h"
#include "H5Cpp.h"

#include "energyGrid.h"
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
   const string inputFile = "../input_file.json";
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

   
   // initialize the enegy grid and eedf
   //
   energyGrid Egrid;
   Egrid.initialize(inputRoot);
   //
   EEDF eedf;
   eedf.initialize(Egrid, inputRoot);   

   
   // load more information and manipulate f0 in time
   //
   vector<vector<double>> vec2d;
   int nrow = 5;
   int ncol = 4;
   vec2d.resize(nrow,vector<double>(ncol,0.0));
   //cout << vec2d.size() << endl;
   //cout << vec2d[0].size() << endl;


   // write energy grid to a specified h5 output file
   // 
   const string outputFile = "output.h5";
   H5File file(outputFile.c_str(), H5F_ACC_TRUNC); // overwrites old
   //
   HDF5dataFile dataFile;
   dataFile.writeOutput(outputFile, "Ecc", Egrid.Ecc, 0);
   dataFile.writeOutput(outputFile, "Ece", Egrid.Ece, 0);
   dataFile.writeOutput(outputFile, "f0", eedf.f0, 1);
   dataFile.writeOutput(outputFile, "Te0", eedf.Te0, 1);
   dataFile.writeOutput(outputFile, "f0", eedf.f0, 1);
   //dataFile.writeOutput(outputFile, "Ecc", Egrid.Ecc, 0); // returns error
   //dataFile.writeOutput(outputFile, "Ecc", Egrid.Ecc, 1); // returns error
   //dataFile.writeOutput(outputFile, "vec2d", vec2d);
   
   cout << "\nEnding simulation" << endl;
   cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" << endl;
   return 1;


   //testing::InitGoogleMock(&argc, argv);
   //return RUN_ALL_TESTS();
}


