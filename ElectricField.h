/***
 * 
 * electric Field class
 *
***/

#ifndef electricField_h
#define electricField_h

using namespace std;

class electricField
{
public:
   double EVpm;   // electric field [V/m]
                        
   void initialize(const Json::Value&, HDF5dataFile&);
};

void electricField::initialize(const Json::Value& root, HDF5dataFile& dataFile)
{
   const Json::Value defValue; // used for default reference
   const Json::Value EField = root.get("ElectricField",defValue);
   if(EField.isObject()) {
      printf("Initializing electric field ...\n");
      Json::Value EVal = EField.get("E",defValue);
      EVpm = EVal.asDouble();
      cout << "E = " << EVpm << " [V/m]" << endl;
   }
   else {
      cout << "value for key \"EField\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   
   dataFile.add(EVpm, "E", 0);
    
   cout << endl;  
}


#endif
