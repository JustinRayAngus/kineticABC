#include "gmock/gmock.h"
#include <iostream>
#include <string>
#include <fstream>
#include "json/json.h"

using namespace std;

int main(int argc, char** argv) {

   // First task is to be able to read energy grid information from
   // a specified input file.

   Json::Value root; // will contain root value after parsing
   Json::Reader reader;
   const Json::Value defValue;  // used for default reference
   string thisInputFile = "input_file.json";
   ifstream ifile(thisInputFile); 
   bool isJsonOK = (ifile !=NULL && reader.parse(ifile, root) );
   if(isJsonOK){
      const Json::Value Egrid = root.get("Egrid",defValue);
      if(Egrid.isObject()){
         printf("\n Reading energy grid information from file: %s %s", 
                    thisInputFile.c_str(), "\n");
         Json::Value EmaxVal = Egrid.get("Emax","");
         double Emax = EmaxVal.asDouble(); // maximum value on energy-space grid
         Json::Value nEVal = Egrid.get("nE","");
         int nE = nEVal.asInt();           // number of energy-space grid points
         double dE = Emax/nE;              // energy-space grid spacing
 
         //cout << "nEVal : " << nEVal.asString() << endl;
         cout << "nE : " << nE << endl;
         cout << "Emax : " << Emax << endl;
         cout << "dE : " << dE << endl;
      }
      else{
         cout << "value for key \"Egrid\" is not object type !" << endl;
      }      

   }
   else 
      cout << "json not OK !!" << endl;

  
   cout << "\n next step is to write stuff to hdf5 file" << endl;
   return 1;


   //testing::InitGoogleMock(&argc, argv);
   //return RUN_ALL_TESTS();
}

