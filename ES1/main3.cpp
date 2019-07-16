
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>		//std::pow();
#include <iomanip>		//std::setw();
#include "../PRNG/BASES/random.h"

using namespace std;

void Blocks(int NoTot, int NoBlk, double *A_i,char* outfile)
{
//   NoTot=1E4;
//   NoBlk=100;
   int TrBlk=int(NoTot/NoBlk);		//Throws per block

   ofstream out(outfile);
   double sigma;
   double true_sum=0, true_sum2=0;
   double true_ave, true_ave2;
   for(int i=0; i<NoBlk; i++)
   {
     true_sum+=A_i[i];						//Vero valor medio					
     true_sum2+=A_i[i]*A_i[i];
     true_ave=true_sum/(double)(i+1);
     true_ave2=true_sum2/(double)(i+1);
     if(i!=0)
       sigma=sqrt((double)(true_ave2-pow(true_ave,2))/i);
     else
       sigma=0;
     out<<i<<setw(15)<<true_ave<<setw(15)<<sigma<<endl;
   }
   out.close();
   return;
}

 
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
   int NoTot=1E6;
   int NoBlk=100;
   int TrBlk=int(NoTot/NoBlk);		//Throws per block

   double L=1.;
   double d=1.5;
   int cont;

   double A_i[NoBlk];
//----------------FIGURA 1--------------------------------

   char outfile1[20]="Es01.3.res";
   for(int i=0; i<NoBlk; i++) 
   {
     cont=0;
     for(int k=0; k<TrBlk; k++)
     {
       double x0=rnd.Rannyu(0., d);
       double theta, x, y;
       do			//Sampling theta uniformly between 0 and pi without using pi's value
       {
         x=rnd.Rannyu(-1., 1.);
         y=rnd.Rannyu(0., 1.);
         theta=acos(x/pow(x*x+y*y,0.5));
       } while (x*x+y*y>=1.);
       if(d<=x0+L*sin(theta)) cont++;		//Checking if bar is onto the line or not
     }
     A_i[i]=2.*L*(TrBlk)/(d*cont);
   }

   Blocks(NoTot, NoBlk, A_i, outfile1);

   rnd.SaveSeed();
   return 0;
}
