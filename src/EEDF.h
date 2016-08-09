/***
 * 
 * electron-energy-distribution function class
 *
***/

#ifndef EEDF_h
#define EEDF_h
//#ifndef __EEDF_H_INCLUDED__
//#define __EEDF_H_INCLUDED__


#include "Gas.h"
#include "energyGrid.h"


using namespace std;

class EEDF
{
public:
   double zeroMom, Te;       // zero moment
   double nunet=0.0;            // dne/dt/ne at n+1/2;
   string type0;         // initial type of EEDF
   vector<double> F0, F0old, F0half; // EEDF
   vector<double> W, D, Flux;  // Energy space adv, diff, and flux at cell-edge
   vector<double> ExcS, IznS; // electronic-excitation source term at cell-center
   
   void initialize(const energyGrid&, const Json::Value&);
   void computeFlux(const Gas&, const energyGrid&, const double&);
   void computeExcS(const Gas&, const energyGrid&);
   void computeIznS(const Gas&, const energyGrid&);
   void advanceF0(const energyGrid&, const double&);

   //double getErDoneTe() const {return Te;}   

private:
   //double Te;
   void setTriDiagCoeffs(vector<double>&, vector<double>&, vector<double>&, vector<double>&,
                         const energyGrid&, const double&);
   void getOpalH(double&, const double&, const double&, 
                 const double&, const double&, const double&);
   void IntExp(double&, const double&, const double&, const double&, const int&);
};


void EEDF::initialize(const energyGrid& Egrid, const Json::Value& root)
{
   F0.assign(Egrid.nE,0.0);
   Flux.assign(Egrid.nE+1,0.0);
   ExcS.assign(Egrid.nE,0.0);
   IznS.assign(Egrid.nE,0.0);
   W.assign(Egrid.nE+1,0.0);
   D.assign(Egrid.nE+1,0.0);
   const Json::Value defValue; // used for default reference
   const Json::Value EEDF = root.get("EEDF",defValue);
   if(EEDF.isObject()) {
      printf("Initializing EEDF ...\n");
      Json::Value TeVal = EEDF.get("Te",defValue);
      Json::Value type = EEDF.get("type",defValue);
      if(TeVal == defValue || type == defValue) {
         printf("ERROR: Te or type not declared in input file\n");
         exit (EXIT_FAILURE);
      } 
      Te = TeVal.asDouble();
      if(Te < 0.0) {
         printf("ERROR: Te is not a positive value in input file\n");
         exit (EXIT_FAILURE);
      }
      type0 = type.asString();
      if(type0=="maxwellian" || type0=="Maxwellian") {
         cout << "Initial EEDF is Maxwellian with Te = " << Te << endl;
         zeroMom = 0;
         for (auto n=0; n<Egrid.nE; n++) {
            F0[n] = 2.0/sqrt(Te*3.14159)/Te*exp(-Egrid.Ecc[n]/Te);
            zeroMom = zeroMom + sqrt(Egrid.Ecc[n])*F0[n]*Egrid.dE;
         }
         transform(F0.begin(), F0.end(), F0.begin(), bind1st(multiplies<double>(),1.0/zeroMom)); 
         //for (auto n=0; n<Egrid.nE; n++) {
         //   F0[n] /= zeroMom;
         //}
         zeroMom = 1.0;
      }
      else {
         cout << "Initial EEDF type = " << type0 << " is not valid " << endl;
         exit (EXIT_FAILURE);
      }
   }
   else {
      cout << "value for key \"Egrid\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   F0old = F0;
   F0half = F0;
  
   //computeExcS(gas, Egrid);
   //computeIznS(gas, Egrid);
   cout << endl;  
}

void EEDF::computeFlux(const Gas& gas, const energyGrid& Egrid, const double& EVpm)
{
   const double Ng = gas.Ng;
   const double mM = gas.mM;
   const double econst = 1.6022e-19;
   const double meconst = 9.1094e-31;
   const double kBconst = 1.3807e-23;
   const double gamma = sqrt(2.0*econst/meconst);
   const vector<double> Ece = Egrid.Ece;
   const int nE = Egrid.nE; // number of cell-center points
   
   // compute advection and diffusion coefficients
   //
   vector<double> PecNum; 
   PecNum.assign(nE+1,0.0);     // Peclet number, |W*dE/D| < 2, 
   
   for (auto n=1; n<nE+1; n++) { // dont need values at E=0
      W[n] = -gamma*Ng*2.0*mM*pow(Egrid.Ece[n],2)*gas.Qelm[n];
      //D[n] = gamma*EVpm*EVpm/3.0/Ng/gas.Qmom[n]*Egrid.Ece[n]
      D[n] = gamma*EVpm*EVpm/3.0*pow(Ece[n],1.5)/(sqrt(Ece[n])*gas.Ng*gas.Qmom[n]+nunet/gamma)
           + gamma*kBconst*Ng*gas.Tg/econst*pow(Egrid.Ece[n],2)*2.0*mM*gas.Qelm[n];
      PecNum[n] = W[n]*Egrid.dE/D[n];
   }
   
   // make sure grid is refined enough for results to be physical
   //
   double PecMax = *max_element(begin(PecNum), end(PecNum));
   double PecMin = *min_element(begin(PecNum), end(PecNum));
   if(PecMax >= 2 || abs(PecMin) >=2) {
      cout << "Maximum Peclet number = " << max(PecMax,abs(PecMin)) << " >=2 " << endl;
      cout << "Grid is not refined enough for results to be physical !!!" << endl;
      exit(EXIT_FAILURE);
   }
   
   // compute Flux = W*F0 - D*dF0/dE
   // boundary conditions are Flux=0 at E=0 and E=Emax
   // E=0 BC is physical and E=Emax BC is to conserve zero moment
   // F0 and dF0/dE are finite at E=0 => W = D = 0 at E=0
   // W*F0=D*DF0/dE at E=Emax
   //
   Flux[0] = 0;
   for (auto n=1; n<nE; n++) { 
      Flux[n] = W[n]*(F0half[n]+F0half[n-1])/2.0 
              - D[n]*(F0half[n]-F0half[n-1])/Egrid.dE;
   }
   Flux[nE] = 0; 
}

void EEDF::computeExcS(const Gas& gas, const energyGrid& Egrid)
{
   const double Ng = gas.Ng;
   const vector<vector<double>> Qexc = gas.Qexc;
   const vector<double> Uexc = gas.Uexc;
   const double econst = 1.6022e-19;
   const double meconst = 9.1094e-31;
   const double gamma = sqrt(2.0*econst/meconst);
   const int nE = Egrid.nE; // number of cell-center points
   
   // loop through each reaction
   //
   double thisU, PL, PU, Esub, deltaE;
   double a, b;
   // deltaQ;
   //double c, d, g, exp1, exp2;
   int thism;
   vector<double> thisQ;
   int numReacs = Uexc.size();
   fill(ExcS.begin(),ExcS.end(),0.0); // reset to zeros;
   for (auto n=0; n<numReacs; n++) {
      thisU = Uexc[n];
      thisQ = Qexc[n];
      //cout << thisQ.size() << endl;
      thism = 1;
      for (auto j=0; j<nE; j++) {
         if( Egrid.Ece[j] < thisU && thisU <= Egrid.Ece[j+1] ) {
         
            // this cell contains thisU
            //
            a = thisU;
            b = Egrid.Ece[j+1];
            if(a==b) a = Egrid.Ece[j]; // prevent divide by zero below using method 2 and 3
 //           deltaQ = thisQ[j+1]-thisQ[j]; // method 2
            PU = Ng*gamma*F0half[j]*(thisQ[j]+thisQ[j+1])/2.0*(b*b-a*a)/2.0; // method 1
 //         PU = Ng*gamma*F0half[j]*(b*b/2.0*(thisQ[j]-deltaQ*a/(b-a))+b*b*b*deltaQ/3.0/(b-a)
 //            -                 a*a/2.0*(thisQ[j]-deltaQ*a/(b-a))-a*a*a*deltaQ/3.0/(b-a));

//            g = log(F0half[j-1]/F0half[j+1])/(Egrid.Ecc[j+1]-Egrid.Ecc[j-1]);           
//            c = thisQ[j]-a*deltaQ/(b-a);
//            d = deltaQ/(b-a); 
//            IntExp(exp1,a,b,g,1);
//            IntExp(exp2,a,b,g,2);          
//            PU = Ng*gamma*F0half[j]*exp(Egrid.Ecc[j]*g)*(c*exp1+d*exp2); // method 3
          
            // Update RHS dint(f)/dE/dt due to source term
            //  
            ExcS[j] -= PU; //Egrid.dE;
            ExcS[0] += PU; //Egrid.dE;
            deltaE = thisU - Egrid.Ece[j];
            //cout << "this J = " << j << endl;
         }
         if(Egrid.Ece[j] >= thisU) {
            //cout << "this j = " << j << endl;
            Esub = Egrid.Ece[j]+deltaE;
  //          deltaQ = thisQ[j+1]-thisQ[j];
            
            // calculate source term for cell j
            //
            a = Egrid.Ece[j];
            b = Esub;
            PL = Ng*gamma*F0half[j]*(thisQ[j]+thisQ[j+1])/2.0*(b*b-a*a)/2.0; // method 1
           // PL = Ng*gamma*F0half[j]*(b*b/2.0*(thisQ[j]-deltaQ*a/Egrid.dE)+b*b*b*deltaQ/3.0/Egrid.dE
           //    -                 a*a/2.0*(thisQ[j]-deltaQ*a/Egrid.dE)-a*a*a*deltaQ/3.0/Egrid.dE);
            a = Esub;
            b = Egrid.Ece[j+1];
            PU = Ng*gamma*F0half[j]*(thisQ[j]+thisQ[j+1])/2.0*(b*b-a*a)/2.0; // method 1
           // PU = Ng*gamma*F0half[j]*(b*b/2.0*(thisQ[j]-deltaQ*a/Egrid.dE)+b*b*b*deltaQ/3.0/Egrid.dE
           //    -                 a*a/2.0*(thisQ[j]-deltaQ*a/Egrid.dE)-a*a*a*deltaQ/3.0/Egrid.dE);

            // begin method 3 (See Hagelaar Eqs 46-51 (note sign error in 51))
            //
//             g = log(F0half[j-1]/F0half[j+1])/(Egrid.Ecc[j+1]-Egrid.Ecc[j-1]);
//             if(j==nE-1) g = log(F0half[j-1]/F0half[j])/(Egrid.Ecc[j+1]-Egrid.Ecc[j]);
//             a = Egrid.Ece[j];
//             b = Esub;
//             c = thisQ[j]-a*deltaQ/Egrid.dE;
//             d = deltaQ/Egrid.dE; 
//             IntExp(exp1,a,b,g,1);
//             IntExp(exp2,a,b,g,2);          
//             PL = Ng*gamma*F0half[j]*exp(Egrid.Ecc[j]*g)*(c*exp1+d*exp2);
//             
//             a = Esub;
//             b = Egrid.Ece[j+1];
//             c = thisQ[j]-a*deltaQ/Egrid.dE;
//             d = deltaQ/Egrid.dE; 
//             IntExp(exp1,a,b,g,1);
//             IntExp(exp2,a,b,g,2); 
//             PU = Ng*gamma*F0half[j]*exp(Egrid.Ecc[j]*g)*(c*exp1+d*exp2);
            //
            // end method 3

            // Update RHS dint(f)/dE/dt due to excitation source term
            //  
            ExcS[j] -= (PL+PU); 
            ExcS[thism-1] += PL; 
            ExcS[thism] += PU;
            thism = thism+1;
         } 
      }
   }
   for (auto j=0; j<nE; j++) {
      ExcS[j] /= Egrid.dE*sqrt(Egrid.Ecc[j]);
   }

}

void EEDF::computeIznS(const Gas& gas, const energyGrid& Egrid)
{
   const double Ng = gas.Ng;
   const vector<vector<double>> Qizn = gas.Qizn;
   const vector<double> Uizn = gas.Uizn;
   const vector<double> energySharing = gas.Sharing;
   const double econst = 1.6022e-19;
   const double meconst = 9.1094e-31;
   const double gamma = sqrt(2.0*econst/meconst);
   const int nE = Egrid.nE; // number of cell-center points
   
   // loop through each ionization reaction
   //
   double thisU, PL, PU, Esub, deltaE;
   double a, b, thisbi0, thiski;
   double thisH, thisNuSink, thisNuSrcP, thisNuSrcS;
   string sharing;
   bool SharedCell = true; // used for equal energy sharing
   int thism, kU, kEU, k2EU;
   vector<double> thisQ, thisIznSrcP, thisIznSrcS;
   int numReacs = Uizn.size();
   fill(IznS.begin(),IznS.end(),0.0); // reset to zeros
   for (auto n=0; n<numReacs; n++) {
      thisU = Uizn[n];
      thisQ = Qizn[n];
      if(energySharing[n]==0) sharing = "zero";
      if(energySharing[n]==1) sharing = "equal";
      if(energySharing[n]==2) {
         sharing = "Opal";
         thisbi0 = gas.bi0[n];
         thiski  = gas.ki[n];
         thisIznSrcP.assign(nE,0.0);
         thisIznSrcS.assign(nE,0.0);
         thisNuSink = 0.0;  // numerical izn rate from sink term
         thisNuSrcP = 0.0;  // numerical izn rate from primary source term
         thisNuSrcS = 0.0;  // numerical izn rate from secondary Source term
      }
      if(energySharing[n]>=3) {
         cout << "ionization energy sharing method must be 0 (zero), 1 (equal), or 2 (Opal)" << endl;
         exit(EXIT_FAILURE);
      }
      //cout << thisQ.size() << endl;
      thism = 1;
      for (auto j=0; j<nE; j++) {
         if( Egrid.Ece[j] < thisU && thisU <= Egrid.Ece[j+1] ) {
         
            // this cell contains thisU
            //
            a = thisU;
            b = Egrid.Ece[j+1];
            PU = Ng*gamma*F0half[j]*(thisQ[j]+thisQ[j+1])/2.0*(b*b-a*a)/2.0;
    
            // Update ionization sink term
            //  
            IznS[j] -= PU;
          
            if(sharing=="equal" || sharing=="zero") {
               IznS[0] += PU; // primaries go to first energy bin
               IznS[0] += PU; // secondaries go to first energy bin
               deltaE = thisU - Egrid.Ece[j];
               //cout << "this J = " << j << endl; 
            }
            else { // source term for Opal method added after sink term loop
               thisNuSink = PU;
               kU = j;      // index for cell that contains thisU
               kEU = kU;    // index for cell containing U+Ecc[0]
               if(thisU+Egrid.Ecc[0]>Egrid.Ece[j+1]) kEU=kU+1;
               k2EU = kU+1; // index for cell containing E+2Ecc[0]
            }
            
         }
         if(Egrid.Ece[j] >= thisU) {
            //cout << "this j = " << j << endl;
            Esub = Egrid.Ece[j]+deltaE;
            
            // calculate sink/source term for cell j
            //
            a = Egrid.Ece[j];
            b = Esub;
            PL = Ng*gamma*F0half[j]*(thisQ[j]+thisQ[j+1])/2.0*(b*b-a*a)/2.0; 
            a = Esub;
            b = Egrid.Ece[j+1];
            PU = Ng*gamma*F0half[j]*(thisQ[j]+thisQ[j+1])/2.0*(b*b-a*a)/2.0;

            // Update ionization sink/source terms for "zero" and "equal"
            //  
            IznS[j] -= (PL+PU); // add local sink term
            if(sharing=="zero") {
               IznS[thism-1] += PL; 
               IznS[thism] += PU;
               IznS[0] += PL+PU; // new electrons go to zero energy
               thism = thism+1;
            }
            if(sharing=="equal") {
               if(SharedCell) {
                  IznS[thism-1] += 2.0*(PL+PU); 
                  SharedCell = false;
               }
               else{
                  IznS[thism-1] += 2.0*PL; 
                  IznS[thism] += 2.0*PU;
                  thism = thism+1;
                  SharedCell = true;
               }
            }
            if(sharing=="Opal") {
               thisNuSink +=(PL+PU);
            }
         } 
      }
      if(sharing=="Opal") { // now update source term for Opal
         for (auto j=0; j<nE; j++) {
           
            // compute source term for primary at this energy
            //
            a = thisU + Egrid.Ecc[j];
            b = min(thisU + 2.0*Egrid.Ecc[j], Egrid.Ece[kEU+1]);
            PU = Ng*gamma*F0half[kEU]*(thisQ[kEU]+thisQ[kEU+1])/2.0*(b*b-a*a)/2.0; 
            getOpalH(thisH, (b+a)/2.0, (b+a)/2.0-Egrid.Ecc[j]-thisU, thisU, thisbi0, thiski);         
            thisIznSrcP[j] += PU*thisH*Egrid.dE; 
            thisNuSrcP += PU*thisH*Egrid.dE;          
            for (auto n=kEU+1; n<k2EU; n++) {
               a = b;
               b = Egrid.Ece[n+1];
               PU = Ng*gamma*F0half[n]*(thisQ[n]+thisQ[n+1])/2.0*(b*b-a*a)/2.0; 
               getOpalH(thisH, (b+a)/2.0, (b+a)/2.0-Egrid.Ecc[j]-thisU, thisU, thisbi0, thiski);  
               thisIznSrcP[j] += PU*thisH*Egrid.dE;          
               thisNuSrcP += PU*thisH*Egrid.dE;
            
            }
            a = b;
            b = min(thisU + 2.0*Egrid.Ecc[j], Egrid.Ece[nE]);
            PU = Ng*gamma*F0half[k2EU]*(thisQ[k2EU]+thisQ[k2EU+1])/2.0*(b*b-a*a)/2.0; 
            getOpalH(thisH, (b+a)/2.0, (b+a)/2.0-Egrid.Ecc[j]-thisU, thisU, thisbi0, thiski);  
            thisIznSrcP[j] += PU*thisH*Egrid.dE; 
            thisNuSrcP += PU*thisH*Egrid.dE;
            
            // compute source term for secondary at this energy
            //
            a = b;
            b = Egrid.Ece[k2EU+1];
            PU = Ng*gamma*F0half[k2EU]*(thisQ[k2EU]+thisQ[k2EU+1])/2.0*(b*b-a*a)/2.0; 
            getOpalH(thisH, (b+a)/2.0, Egrid.Ecc[j], thisU, thisbi0, thiski);  
            thisIznSrcS[j] += PU*thisH*Egrid.dE; 
            thisNuSrcS += PU*thisH*Egrid.dE;
            for (auto n=k2EU+1; n<nE; n++) {
               a = b;
               b = Egrid.Ece[n+1];
               PU = Ng*gamma*F0half[n]*(thisQ[n]+thisQ[n+1])/2.0*(b*b-a*a)/2.0; 
               getOpalH(thisH, (b+a)/2.0, Egrid.Ecc[j], thisU, thisbi0, thiski);  
               thisIznSrcS[j] += PU*thisH*Egrid.dE;
               thisNuSrcS += PU*thisH*Egrid.dE;         
            }                
            
            // update indices for integration limits
            // 
            kEU  = kEU+1; // lower index for primary source integral
            k2EU = min(k2EU+2,Egrid.nE-1); // upper index for primary source integral
            if(kEU==Egrid.nE) break;
            
         }
          
         // need to renormalize both sources such that integrals are equal -sink
         // then add them to the total
         //
         for (auto j=0; j<nE; j++) {
            IznS[j] += thisNuSink/thisNuSrcP*thisIznSrcP[j];
            IznS[j] += thisNuSink/thisNuSrcS*thisIznSrcS[j];
         }
//       transform(thisIznSrcP.begin(), thisIznSrcP.end(), thisIznSrcP.begin(), 
//                 bind1st(multiplies<double>(),thisNuSink/thisNuSrcP)); 
//       transform(thisIznSrcS.begin(), thisIznSrcS.end(), thisIznSrcS.begin(), 
//                 bind1st(multiplies<double>(),thisNuSink/thisNuSrcS)); 
//       transform(IznS.begin(), IznS.end(), IznS.begin(),
//                 bind1st(adds<double>(),thisIznSrcP)); 
//       transform(IznS.begin(), IznS.end(), IznS.begin(),
//                 bind1st(adds<double>(),thisIznSrcS)); 

      }   
                    
   }
   
   
   nunet = 0.0;
   double zeroMomHalf = 0.0;
   for (auto j=0; j<nE; j++) {
      nunet += IznS[j];
      IznS[j] /= Egrid.dE*sqrt(Egrid.Ecc[j]);
      //nunet += IznS[j]*sqrt(Egrid.Ecc[j])*Egrid.dE;
      zeroMomHalf += F0half[j]*sqrt(Egrid.Ecc[j])*Egrid.dE; // should be unity
   }
   nunet /=zeroMomHalf; // for numerical conservation purposes
}

void EEDF::advanceF0(const energyGrid& Egrid, const double& dt)
{
   const int nE = Egrid.nE;
   /*
   // Explicit forward advance
   //
   //for (auto n=0; n<nE; n++) {
   //   F0[n] = F0[n] - dt*(Flux[n+1]-Flux[n])/Egrid.dE/sqrt(Egrid.Ecc[n])
   //                 + dt*(ExcS+IznS);
   //}
   */
   
   // Implicit Crank-Nicolson advance
   //
   
   // Compute the tridiagonal matrix components
   //
   vector<double> a(nE,0.0), b(nE,0.0), c(nE,0.0), d(nE,0.0);
   setTriDiagCoeffs(a,b,c,d,Egrid,dt);
   
   // add inelastic source/sink term to d vector
   //
   for (auto n=0; n<nE; n++) {
      d[n] += 4.0*dt*(ExcS[n]+IznS[n]-nunet*F0half[n]);
   }
   
   // Row reduce
   //
   for (auto n=1; n<nE; n++) {
      b[n] = b[n] - a[n]*c[n-1]/b[n-1];
      d[n] = d[n] - a[n]*d[n-1]/b[n-1];
   }
   
   // Back solve
   //
   F0[nE-1] = d[nE-1]/b[nE-1];
   for (auto n=nE-2; n>-1; n--) {
      F0[n] = (d[n] - F0[n+1]*c[n])/b[n]; 
      if(F0[n]<0) {
         cout << "F0("<< n <<")= " << F0[n] << endl;
         //F0[n] = 1.0e-15;
      }
   }
   
   // update moments and F0half
   //
   Te = 0;
   zeroMom = 0;
   for (auto n=0; n<nE; n++) {
      F0half[n] = (F0[n]+F0old[n])/2.0;
      Te = Te + 2.0/3.0*pow(Egrid.Ecc[n],1.5)*F0[n]*Egrid.dE;
      zeroMom = zeroMom + sqrt(Egrid.Ecc[n])*F0[n]*Egrid.dE;
   }
   // For numerical conservation purposes divide by numerical zeroMom, which should be close to unity 
   // Finite volume method should not need this...I need to write out the matrix and look at it
   //
   //transform(F0.begin(), F0.end(), F0.begin(), bind1st(multiplies<double>(),1.0/zeroMom)); 

}

void EEDF::IntExp(double& soln, const double& a, const double& b, const double& g, const int& pow)
{ 
   if(pow==1) {
      soln = -( exp(-b*g)*(b*g+1.0) - exp(-a*g)*(a*g+1.0) )/(g*g);
   }
   if(pow==2) {
      soln = -( exp(-b*g)*(b*b*g*g + 2.0*(b*g+1.0) ) 
           -    exp(-a*g)*(a*a*g*g + 2.0*(a*g+1.0) ) )/(g*g*g);
   }
   if(pow!=1 && pow!=2) {
      cout << "Error: call to IntExp only accepts pow=1 and 2" << endl;
      exit(EXIT_FAILURE);
   }

}

void EEDF::getOpalH(double& H, const double& PrimaryE, const double& SecondaryE, 
                    const double& Ui, const double& bi0, const double& ki)
{ 
   
   // H normalized such that int_{U}^{(PrimaryE-U)/2}H(PrimaryE,SecondaryE)dSecondaryE = 1   
   double bi = bi0;
   if(PrimaryE>=exp(ki)) bi = bi0*ki/log(PrimaryE);
   if(PrimaryE<=Ui) {
      cout << " PrimaryE <= Ui in call to Opal H !!! " << endl;
      exit(EXIT_FAILURE);
   }   
   H = bi/atan((PrimaryE-Ui)/(2.0*bi))/(bi*bi + SecondaryE*SecondaryE);

}

void EEDF::setTriDiagCoeffs(vector<double>& a, vector<double>& b, vector<double>& c, 
                            vector<double>& d, const energyGrid& Egrid, const double& dt)
{  
   // create vectors of Courant and alphas
   //
   const int nE = Egrid.nE;
   vector<double> couraM(nE,0.0), alphaM(nE,0.0); // W and D at n-1/2
   vector<double> couraP(nE,0.0), alphaP(nE,0.0); // W and D at n+1/2
      
   for(auto n=0; n<nE; n++) {
      couraM[n] = W[n]*dt/Egrid.dE/sqrt(Egrid.Ecc[n]);
      alphaM[n] = D[n]*dt/(Egrid.dE*Egrid.dE*sqrt(Egrid.Ecc[n]));
      couraP[n] = W[n+1]*dt/Egrid.dE/sqrt(Egrid.Ecc[n]);
      alphaP[n] = D[n+1]*dt/(Egrid.dE*Egrid.dE*sqrt(Egrid.Ecc[n]));  
   }
   
   // set tridiag coefficients 
   //
   for(auto n=0; n<nE; n++) {
      a[n] = - (couraM[n] + 2.0*alphaM[n]);
      b[n] = 4.0 + couraP[n] - couraM[n] + 2.0*(alphaP[n] + alphaM[n]);
      c[n] = couraP[n] - 2.0*alphaP[n];     
   }
   //const double EmaxBC = 0; // zero value at Emax+1/2
   const double EmaxBC = (2.0*D[nE]+Egrid.dE*W[nE])/(2.0*D[nE]-Egrid.dE*W[nE]); // 0 Flux at Emax
   b[nE-1] += c[nE-1]*EmaxBC;
      
   // set d coefficient A*F0^{n+1} = d(F0^{n})
   //
   d[0] = (4.0 - couraP[0] - 2.0*alphaP[0])*F0old[0] 
        - (couraP[0] - 2.0*alphaP[0])*F0old[1];
   for(auto n=1; n<nE-1; n++) {
      d[n] = (couraM[n] + 2*alphaM[n])*F0old[n-1]
           + (4.0 + couraM[n] - couraP[n] - 2.0*(alphaM[n] + alphaP[n]))*F0old[n]
           - (couraP[n] - 2.0*alphaP[n])*F0old[n+1];
   }
   d[nE-1] = (couraM[nE-1]+2.0*alphaM[nE-1])*F0old[nE-2]
           + (4.0+couraM[nE-1]-couraP[nE-1]-2.0*(alphaM[nE-1]+alphaP[nE-1]))*F0old[nE-1];
   d[nE-1] -= (couraP[nE-1]-2.0*alphaP[nE-1])*EmaxBC*F0old[nE-1];
}

#endif
