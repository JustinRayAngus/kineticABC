#include "gmock/gmock.h"
#include <iostream>

using namespace std;

int main(int argc, char** argv) {

   int i;
   cout << "Please enter an integer value: ";
   cin >> i;
   cout << "The value you entered is " << i << ".\n";
   cout << "Next is learn how to read input files using json \n";
   cout << "and write output files to hdf5. \n";

   return 0;


   //testing::InitGoogleMock(&argc, argv);
   //return RUN_ALL_TESTS();
}

