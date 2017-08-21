/***
 * 
 * background Gas class
 *
***/

#ifndef Gas_h
#define Gas_h
//#ifndef __Gas_H_INCLUDED__
//#define __Gas_H_INCLUDED__

#include "energyGrid.h"


using namespace std;

class Gas
{
public:
   double Tg, Pg, Ng;   // gas temp [K], press [Torr], and dens [1/m^3]
   string gasName;      // gas species name
   vector<double> Qelm; // elastic-momentum-transfer cross section [m^2]
   double mM;           // electron mass / gas species mass
   vector<vector<double>> Qexc; // electronic-excitation cross section [m^2]
   vector<double> Uexc;         // transition thresholds [eV]
   //vector<double> Qvib; // vibrational-excitation cross section [m^2]
   //double Uvib;
   vector<vector<double>> Qizn; // ionization cross section [m^2]
   vector<double> Uizn, Sharing, bi0, ki;
   vector<double> Qmom; // momentum transfer cross section
               
                        
   void initialize(const energyGrid&, const Json::Value&);

private:
   string xsecsFile;     // xsecs file
   
   void loadXsecs(const energyGrid&, const string&);
   void interpXsecs(const vector<double>&, vector<double>&,
                    const vector<double>&, const vector<double>&,
                    const int&, const double&);
};

void Gas::initialize(const energyGrid& Egrid, const Json::Value& root)
{
   Qelm.assign(Egrid.nE+1,0.0);
   Qexc.resize(1,vector<double>(Egrid.nE+1,0.0));
   Uexc.assign(1,0.0);
   Qizn.resize(1,vector<double>(Egrid.nE+1,0.0));
   Uizn.assign(1,0.0);
   Sharing.assign(1,0.0);
   bi0.assign(1,0.0);
   ki.assign(1,0.0);
   Qmom.assign(Egrid.nE+1,0.0);
   const Json::Value defValue; // used for default reference
   const Json::Value Gas = root.get("Gas",defValue);
   if(Gas.isObject()) {
      printf("Initializing background gas ...\n");
      Json::Value TgVal = Gas.get("Tg",defValue);
      Json::Value PgVal = Gas.get("Pg",defValue);
      Json::Value nameVal = Gas.get("name",defValue);
      Json::Value xsecsFileVal = Gas.get("xsecsFile",defValue);
      if(TgVal == defValue || PgVal == defValue || 
         nameVal == defValue || xsecsFileVal == defValue) {
         printf("ERROR: 'Gas' section not specified correctly in input file\n");
         exit (EXIT_FAILURE);
      } 
      Tg = TgVal.asDouble(); // [K]
      Pg = PgVal.asDouble();
      Ng = Pg/760*273/Tg*2.6868e25; // gas density [1/m^3]
      if(Tg < 0.0) {
         printf("ERROR: Tg is not a positive value in input file\n");
         exit (EXIT_FAILURE);
      }
      gasName = nameVal.asString();
      xsecsFile = xsecsFileVal.asString();
      if(gasName=="nitrogen" || gasName=="Nitrogen" || 
         gasName=="helium" || gasName=="Helium" ) {
         
         double atomicM;
         if(gasName=="nitrogen" || gasName=="Nitrogen") {
             atomicM = 2.0*14.007; // atomic mass of nitrogen molecule
         }
         if(gasName=="helium" || gasName=="Helium") {
             atomicM = 4.0062;     // atomic mass of helium atom
         }

         const double meconst = 9.1094e-28;   // electron mass [g]
         const double amuconst = 1.6605e-24;  // atomic mass unit [g]
         mM = meconst/(amuconst*atomicM); // me/M
         
         cout << "Background gas is " << gasName << " ... " << endl;
         cout << "Tg = " << Tg << " K" << endl;
         cout << "Pg = " << Pg << " Torr" << endl;
         cout << "Ng = " << Ng << " 1/m^3" << endl;
         cout << "m/M = " << mM << endl;
         cout << "xsecs file = " << xsecsFile << endl;
         loadXsecs(Egrid, xsecsFile);
         
         // add all xsecs together to get Qmom
         //
         const int nQ = Qmom.size();
         const int nExc = Qexc.size();
         vector<double> thisQexc;
         Qmom = Qelm;
         for (auto m=0; m<nExc; m++) {
            thisQexc = Qexc[m];
            for (auto n=0; n<nQ; n++) {
               Qmom[n] += thisQexc[n];
            }
         }
         const int nIzn = Qizn.size();
         vector<double> thisQizn;
         for (auto m=0; m<nIzn; m++) {
            thisQizn = Qizn[m];
            for (auto n=0; n<nQ; n++) {
               Qmom[n] += thisQizn[n];
            }
         }     
         
      }
      else {
         cout << "Gas species " << gasName << " is not a valid option" << endl;
         exit (EXIT_FAILURE);
      }
   }
   else {
      cout << "value for key \"Gas\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   
   cout << endl;  
}

void Gas::loadXsecs(const energyGrid& Egrid, const string& xsecsFile)
{
   cout << "Loading cross section data ... " << endl;
   ifstream xfile(xsecsFile.c_str());
   int firstExc = 1;
   int firstIzn = 1;
   if(xfile.is_open()) {
      double mM2, thisE, thisQ;
      double thisU, thisTransition, createInvQ;
      double thisSharing, thisbi0, thiski;
      vector<double> Etemp(1,0.0), Qtemp(1,0.0);
      string str;
      size_t pos, posExc, posIzn;
      while (getline(xfile,str)) {
      
         // look for elastic reactions
         //
         pos = str.find("ELASTIC");
         if(pos!=string::npos) {
           // cout << str << endl;
            getline(xfile,str);
            cout << "Elastic reaction: " << str << endl;
            xfile >> mM2;
            //cout << "m/M = " << mM2 << endl; 
            getline(xfile,str); // parse m/M line
            getline(xfile,str); // parse 1. 1. line
            getline(xfile,str); // parse ----- line
            //int count = 0;
            while (xfile >> thisE && xfile >> thisQ) {
               if(thisE==0) {
                  Qtemp[0] = thisQ;
               } 
               else {
                  Etemp.push_back(thisE);
                  Qtemp.push_back(thisQ);
               }
            }
            xfile.clear();
            xfile.ignore('\n');
            //cout << Etemp.back() << endl;
            //cout << "Qtemp size = " << Qtemp.size() << endl;
            //cout << "Etemp[0] = " << Etemp[0] << endl;
            //cout << "Etemp[1] = " << Etemp[1] << endl;
            //cout << "Etemp[2] = " << Etemp[2] << endl;
            interpXsecs(Egrid.Ece, Qelm, Etemp, Qtemp, 1, 0);
            Etemp.clear();
            Qtemp.clear();
            //Qmom = Qelm;
         }
         
         // look for Excitation reactions
         //
         posExc = str.find("EXCITATION");
         if(posExc!=string::npos) {
            //cout << str << endl;
            getline(xfile,str);
            cout << "Excitation reaction: " << str << endl;;
            xfile >> thisU;
            //cout << "energy threshold: " << thisU << endl;
            xfile >> thisTransition;
            if(thisTransition==1) {
               //cout << "Transition is allowed" << endl;
            }
            if(thisTransition==0) {
               //cout << "Transition is forbidden" << endl;
            }
            if(thisTransition!=0 && thisTransition!=1) {
               cout << "Transition type not specified, must be 0 or 1" << endl;
               exit (EXIT_FAILURE);
            }
            xfile >> createInvQ;
            if(createInvQ==1) {
               //cout << "Inverse xsec will be created" << endl;
            }
            if(createInvQ==0) {
               //cout << "Inverse xsec will not be created" << endl;
            }
            if(createInvQ!=0 && createInvQ!=1) {
               cout << "Create cross section flag not set, must be 0 or 1" << endl;
               exit (EXIT_FAILURE);
            }
            getline(xfile,str);
            getline(xfile,str);
            getline(xfile,str);
            //cout << str << endl;
            while (xfile >> thisE && xfile >> thisQ) {
               Etemp.push_back(thisE);
               Qtemp.push_back(thisQ);
            }
            //cout << Qtemp.front() << endl;
            //cout << Qtemp.back() << endl;
            xfile.clear();
            xfile.ignore('\n');
            vector<double> thisQexc(Egrid.nE+1,0.0);
            interpXsecs(Egrid.Ece, thisQexc, Etemp, Qtemp, thisTransition, thisU);
            Etemp.clear();
            Qtemp.clear();
            if(firstExc) {
               Qexc[0] = thisQexc;
               Uexc[0] = thisU;
               firstExc = 0;
            }
            else {
               Qexc.push_back(thisQexc);
               Uexc.push_back(thisU);
            }
         }
         
         // look for Ionization reactions
         //
         posIzn = str.find("IONIZATION");
         if(posIzn!=string::npos) {
            //cout << str << endl;
            getline(xfile,str);
            cout << "Ionization reaction: " << str << endl;;
            xfile >> thisU;
            //cout << "energy threshold: " << thisU << endl;
            xfile >> thisSharing;
            if(thisSharing==0) cout << " zero energy sharing method " << endl;
            if(thisSharing==1) cout << " equal energy sharing method " << endl;
            if(thisSharing==2) cout << " Opal energy sharing method " << endl;
            if(thisSharing!=0 && thisSharing!=1 && thisSharing!=2) {
               cout << " Sharing method must be 0 (zero), 1 (equal), or 2 (Opal) " << endl;
               exit (EXIT_FAILURE);
            }
            xfile >> thisbi0;
            cout << "bi0 = " << thisbi0 << endl;
            xfile >> thiski;
            cout << "ki = " << thiski << endl;
            getline(xfile,str);
            getline(xfile,str);
            getline(xfile,str);
            //cout << str << endl;
            while (xfile >> thisE && xfile >> thisQ) {
               Etemp.push_back(thisE);
               Qtemp.push_back(thisQ);
            }
            //cout << Qtemp.front() << endl;
            //cout << Qtemp.back() << endl;
            xfile.clear();
            xfile.ignore('\n');
            vector<double> thisQizn(Egrid.nE+1,0.0);
            interpXsecs(Egrid.Ece, thisQizn, Etemp, Qtemp, 1, thisU);
            Etemp.clear();
            Qtemp.clear();
            if(firstIzn) {
               Qizn[0] = thisQizn;
               Uizn[0] = thisU;
               Sharing[0] = thisSharing;
               bi0[0] = thisbi0;
               ki[0] = thiski;
               firstIzn = 0;
            }
            else {
               Qizn.push_back(thisQizn);
               Uizn.push_back(thisU);
               Sharing.push_back(thisSharing);
               bi0.push_back(thisbi0);
               ki.push_back(thiski);
            }
         }
         
         
      }
   }
   else {
      cout << "ERROR: could not open xsecs file" << endl;
      exit (EXIT_FAILURE);
   }
         //   cout << Qexc.size() << endl;
         //   cout << Qexc[Qexc.size()-1].size() << endl;
   
}

void Gas::interpXsecs(const vector<double>& Ece, vector<double>& Q, 
                      const vector<double>& Edata, const vector<double>& Qdata,
                      const int& allowed, const double& U)
{
   int thism = 0;
   int nmax = Ece.size();
   double thisE, a, b;
   for (auto n=0; n<nmax; n++) {
      thisE = Ece[n];
      if (thisE <= max(Edata.front(),U) && U==0.0) { // elastic
         Q[n] = Qdata.front()*thisE/(Edata.front()+1e-20);
         if (Edata.front()==0 && thisE==0) {
            Q[n] = Qdata.front();
         }
      }
      if (thisE <= max(Edata.front(),U) && U!=0.0) { // inelastic
         Q[n] = 0.0;
      }
      if (thisE > max(Edata.front(),U) && thisE < Edata.back()) {
         //thism = 0;
         while (Edata[thism] <= thisE) {
           thism = thism+1;
         }
         if (Qdata[thism-1]==0 || Qdata[thism]==0) { // linear interp
            a = (Edata[thism]-thisE)/(Edata[thism]-max(Edata[thism-1],U));
            b = 1.0-a;
            Q[n] = b*Qdata[thism] + a*Qdata[thism-1];
         }
         else { // log interp
            a = log10(Edata[thism]/thisE)/log10(Edata[thism]/(max(Edata[thism-1],U)));
            b = 1.0-a;
            Q[n] = b*log10(Qdata[thism]) + a*log10(Qdata[thism-1]);
            Q[n] = pow(10,Q[n]); // linear interp on log10 scale
         }
      }
      if (thisE >= Edata.back()) {
      //cout << "ARE WE HERE" << endl;
         if (allowed) { // Q ~ ln(E)/E
            Q[n] = Qdata.back()*log(thisE)/log(Edata.back())*Edata.back()/thisE;
         }
         else {  // Q ~ 1/E^3
            Q[n] = Qdata.back()*pow(Edata.back()/thisE,3);
         }
      }         
   } 
}


#endif
