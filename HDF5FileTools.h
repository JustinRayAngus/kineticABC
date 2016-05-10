/***
 *
 * HDF5 File Tools class
 *
***/

#ifndef HDF5FileTools_h
#define HDF5FileTools_h

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "H5Cpp.h"

using namespace std;
#ifndef H5_NO_NAMESPACE
   using namespace H5;
#endif

class HDF5FileTools
{
public:
   void writeOutput(const string&, const string&, const vector<double>&);
};

void HDF5FileTools::writeOutput(const string& outputFile, const string& varName,
                                const vector<double> &varData)
{                             
   // first have to copy the vector data to an array for writing purposes
   //
   const int Nvar = varData.size(); 
   const int RANK = 1;
   //cout << varData.size() << endl;
   double data[Nvar];
   for (int i=0; i<Nvar; i++)
   {
      data[i] = varData[i];
   } 

   // now write data to hdf5 file and look for errors
   //
   try
   {
      Exception::dontPrint();
      H5File file( outputFile.c_str(), H5F_ACC_RDWR );
      hsize_t dimsf[1];
      dimsf[0] = varData.size();
      DataSpace dataspace( RANK, dimsf);
      FloatType datatype (PredType::NATIVE_DOUBLE);
      datatype.setOrder(H5T_ORDER_LE);
      DataSet dataset = file.createDataSet( varName.c_str(), datatype, dataspace);
      dataset.write( data, PredType::NATIVE_DOUBLE); 
      cout << "Variable " << varName << " written to file: " << outputFile << endl; 
   
      dataset.close();
      dataspace.close();
      file.close();
      //delete data;
   } 
   catch( FileIException error )
   {
      error.printError();
      //return -1;
   }
   catch( DataSetIException error )
   {
      error.printError();
      //return -1;
   }
   catch( DataSpaceIException error )
   {
      error.printError();
      //return -1;
   }
   catch( DataTypeIException error )
   {
      error.printError();
      //return -1;
   }

}


#endif


