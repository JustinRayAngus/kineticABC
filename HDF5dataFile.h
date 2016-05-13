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
   void writeOutput(const string&, const char*, const double&, const bool&);
   void writeOutput(const string&, const char*, const vector<double>&, const bool&);
   
private:
   bool varExists(const H5File&, const char*);
};

void HDF5dataFile::writeOutput(const string& outputFile, const char* varName,
                               const double& varData, const bool& growVar)
{      
   // Open file and check to see if output variable already exists
   //
   H5File file(outputFile.c_str(), H5F_ACC_RDWR);
   if(varExists(file,varName) && !growVar) {
      cout << "ERROR: Output variable " << varName << " already exists " << endl;
      cout << "       and dataset is not selected to grow " << endl;
	  exit (EXIT_FAILURE);
   }
                      
   // Write data to hdf5 file and try to catch errors
   //
   //H5S_class_t type = H5S_SCALAR;
   PredType varType = PredType::NATIVE_DOUBLE;
   //Exception::dontPrint();   
   DataSpace dataspace(H5S_SCALAR);
   DataType datatype(varType);
   DataSet dataset = file.createDataSet(varName, datatype, dataspace);
   dataset.write(&varData, varType); 
   cout << "Scalar " << varName << " written to file: " << outputFile << endl; 
   
   // close opened stuff
   //
   dataset.close();      
   dataspace.close();
   file.close();
}


void HDF5dataFile::writeOutput(const string& outputFile, const char* varName,
                               const vector<double>& varData, const bool& growVar)
{     
   // Open file and check to see if output variable already exists
   //
   H5File file(outputFile.c_str(), H5F_ACC_RDWR);
   if(varExists(file,varName) && !growVar) {
      cout << "ERROR: Output variable " << varName << " already exists " << endl;
      cout << "       and dataset is not selected to grow " << endl;
	  exit (EXIT_FAILURE);
   }
   
   // Need to copy vectors to an array for writing purposes
   //
   int RANK;
   PredType varType = PredType::NATIVE_DOUBLE;
   const int Nvar = varData.size(); //cout << varData.size() << endl;
   double data[Nvar];
   for (int i=0; i<Nvar; i++) {
      data[i] = varData[i];
   } 
   
   if(!varExists(file,varName)) { // then create it
      
      if(growVar) { // create extendable dataset
      
         RANK = 2;
         hsize_t dimsf[RANK];
         dimsf[0] = Nvar; 
         dimsf[1] = 1;
         hsize_t mdimsf[RANK] = {varData.size(), H5S_UNLIMITED};
         DataSpace mdataspace(RANK, dimsf, mdimsf); // memory dataspace
         
         DSetCreatPropList cparms;
         hsize_t chunk_dims[2]={varData.size(),1}; // chunk size doesn't effect here
         cparms.setChunk(RANK,chunk_dims);
         
         DataType datatype(varType);
         DataSet dataset = file.createDataSet(varName, datatype, mdataspace, cparms);
         
         hsize_t size[RANK] = {varData.size(), 1};
         dataset.extend(size);
         DataSpace dataspace = dataset.getSpace();
         hsize_t offset[RANK] = {0,0};
         hsize_t dims[RANK] = {varData.size(), 1};
         dataspace.selectHyperslab(H5S_SELECT_SET,dims,offset); // data dataspace
                  
         dataset.write(data, varType, mdataspace, dataspace);
         cout << "Extendable Vector " << varName << " written to file: " 
              << outputFile << endl;    
              
         // close opened stuff
         //
         dataset.close();      
         dataspace.close();
         file.close();   
      }
      else { // create non-extendable dataset
      
         RANK = 1;
         hsize_t dimsf[RANK];
         dimsf[0] = Nvar; //varData.size();
         DataSpace dataspace(RANK, dimsf);
         DataType datatype(varType);
         DataSet dataset = file.createDataSet(varName, datatype, dataspace);
         dataset.write(data, varType);
         cout << "Non-extendable vector " << varName << " written to file: "
              << outputFile << endl;   
              
         // close opened stuff
         //
         dataset.close();      
         dataspace.close();
         file.close();     
      }
   }
   else { // append to extendable dataset
   
      // open extendable dataset and get ID's to stuff
      //
      RANK = 2;
      DataSet dataset = file.openDataSet(varName);
      hid_t fileID, dsetID, dspaceID;
      fileID = H5Fopen(outputFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      dsetID = H5Dopen(fileID, varName, H5P_DEFAULT);
      dspaceID = H5Dget_space(dsetID);
      
      // get dimensions of extendable dataset
      //
      const int ndims = H5Sget_simple_extent_ndims(dspaceID);
      if(ndims == 1) {
         cout << "ERROR: Trying to extend non-extendable variable " 
              << varName << endl;
	     exit (EXIT_FAILURE);
      }
      hsize_t dims[ndims];
      H5Sget_simple_extent_dims(dspaceID,dims,NULL); // sets dims
      
      // extend dataset and write
      //
      hsize_t offset[RANK], size[RANK];
      offset[0] = 0;
      offset[1] = dims[1]; // equal to total number of previous time steps
      hsize_t dimsappend[RANK];
      dimsappend[0] = Nvar;
      dimsappend[1] = 1;
      dims[1] = dims[1] + dimsappend[1];
      size[0] = dims[0];
      size[1] = dims[1];
      dataset.extend(size);
      DataSpace fspace = dataset.getSpace();
      fspace.selectHyperslab(H5S_SELECT_SET, dimsappend, offset);
      DataSpace mspace(RANK,dimsappend);
      dataset.write(data, varType, mspace, fspace);
      cout << "Extendable Vector " << varName << " updated in file: " 
           << outputFile << endl; 
        
      // close opened stuff
      //
      dataset.close();      
      fspace.close();
      mspace.close();
      file.close();
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


