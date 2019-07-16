/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>		//std::pow();
#include <iomanip>		//std::setw();
#include "../PRNG/BASES/random.h"


using namespace std;
 
int main (int argc, char *argv[]){

//----------------STARTING PROCESS----------------
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("../PRNG/BASES/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../PRNG/BASES/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

//----------------NOW YOU CAN MODIFY----------------
   int NoTot=1E4;
   int N[4]={1, 2, 10, 100};
   double x[4][NoTot];
   int Npin=100;

//----------------UNIFORM DICE----------------------------

   double xmin1=1.;
   double xmax1=6.;
   double h=(xmax1-xmin1)/(double)Npin;
   ofstream out("Es01.2.1.res");

   for(int j=0; j<NoTot; j++) {for(int k=0; k<4; k++) x[k][j]=0;}

   for(int k=0; k<4; k++)
   {
     for(int i=0; i<N[k]; i++)
     {
       for(int j=0; j<NoTot; j++)  x[k][j]+=int(rnd.Rannyu(0.,6)+1)/(double)(N[k]);
     }
   }

   for(int j=0; j<NoTot; j++)
   {
     for(int k=0; k<4; k++) out<<x[k][j]<<setw(15);
     out<<endl;
   }
   
   out.close();

//----------------EXPONENTIAL DICE----------------------------

   double xmin2=0.;
   double xmax2=8.;
   h=(xmax1-xmin1)/(double)Npin;
   ofstream out2("Es01.2.2.res");

   for(int j=0; j<NoTot; j++) {for(int k=0; k<4; k++) x[k][j]=0;}

   for(int k=0; k<4; k++)
   {
     for(int i=0; i<N[k]; i++)
     {
       for(int j=0; j<NoTot; j++)  x[k][j]+=rnd.Exp(1.)/(double)(N[k]);
     }
   }

   for(int j=0; j<NoTot; j++)
   {
     for(int k=0; k<4; k++) out2<<x[k][j]<<setw(15);
     out2<<endl;
   }
   out2.close();

//----------------LORENTIAN DICE-------------------------------

   double xmin3=-10.;
   double xmax3=10.;
   h=(xmax1-xmin1)/(double)Npin;
   ofstream out3("Es01.2.3.res");

   for(int j=0; j<NoTot; j++) {for(int k=0; k<4; k++) x[k][j]=0;}

   for(int k=0; k<4; k++)
   {
     for(int i=0; i<N[k]; i++)
     {
       for(int j=0; j<NoTot; j++)  x[k][j]+=rnd.Lorentz(0., 1.)/(double)(N[k]);
     }
   }

   for(int j=0; j<NoTot; j++)
   {
     for(int k=0; k<4; k++) out3<<x[k][j]<<setw(15);
     out3<<endl;
   }
   out3.close();

   rnd.SaveSeed();
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
