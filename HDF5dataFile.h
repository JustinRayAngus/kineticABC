/***
 *
 * HDF5 data file class
 *
***/

#ifndef HDF5dataFile_h
#define HDF5dataFile_h

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "H5Cpp.h"

using namespace std;
#ifndef H5_NO_NAMESPACE
   using namespace H5;
#endif

class HDF5dataFile
{
public:
   void writeOutput(const string&, const char*, const vector<double>&);

private:
   bool varExists(const H5File&, const char*);
};

void HDF5dataFile::writeOutput(const string& outputFile, const char* varName,
                                const vector<double> &varData)
{ 
     
   // Open file and check to see if output variable already exists
   //
   H5File file(outputFile.c_str(), H5F_ACC_RDWR);
   bool growVar = 0;
   if(varExists(file,varName) && !growVar) {
      cout << "ERROR: Output variable " << varName << " already exists " << endl;
      cout << "       and dataset is not selected to grow " << endl;
	  exit (EXIT_FAILURE);
   }
                      
   // Need to copy vectors to an array for writing purposes
   //
   const int Nvar = varData.size(); //cout << varData.size() << endl;
   const int RANK = 1;
   double data[Nvar];
   for (int i=0; i<Nvar; i++) {
      data[i] = varData[i];
   } 

   // Write data to hdf5 file and try to catch errors
   //
   hsize_t dimsf[RANK];
   dimsf[0] = varData.size();
   PredType varType = PredType::NATIVE_DOUBLE;
   try {
      Exception::dontPrint();
      
      DataSpace dataspace(RANK, dimsf);
      DataType datatype(varType);
      DataSet dataset = file.createDataSet(varName, datatype, dataspace);
      dataset.write(data, varType); 
      cout << "Variable " << varName << " written to file: " << outputFile << endl; 
      
      dataset.close();      
      dataspace.close();
      file.close();
      //delete file;
   } 
   catch(FileIException error) {
      error.printError();
      exit (EXIT_FAILURE);
   }
   catch(DataSetIException error) {
      error.printError();
      exit (EXIT_FAILURE);
   }
   catch(DataSpaceIException error) {
      error.printError();
      exit (EXIT_FAILURE);
   }
   catch(DataTypeIException error) {
      error.printError();
      exit (EXIT_FAILURE);
   }

}

bool HDF5dataFile::varExists(const H5File& file2, const char* thisVar)
{
   try {
      Exception::dontPrint();
	  DataSet dataset0 = DataSet(file2.openDataSet(thisVar));
	  dataset0.close();
	  return 1; // dataset does not exist in output file yet
   }
   catch( ... ) {  
      return 0; // dataset already exists in output file
   }
}

#endif


