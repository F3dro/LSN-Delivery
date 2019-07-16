
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>		//std::pow();
#include <iomanip>		//std::setw();
#include "../PRNG/BASES/random.h"

using namespace std;

void Blocks(int NoTot, int NoBlk, double *meanval,char* outfile)
{
//   NoTot=1E4;
//   NoBlk=100;
   int TrBlk=int(NoTot/NoBlk);		//Throws per block

   ofstream out(outfile);
   double ave;
   double sigma;
   double true_sum=0, true_sum2=0;
   double true_ave, true_ave2;
   for(int i=0; i<NoBlk; i++)
   {
     double sum=0;
     for(int j=0; j<TrBlk; j++)
        sum+=meanval[i*NoBlk+j];		//Integral to be computed
     ave=(double)sum/(double)TrBlk;		//Media su ogni blocco (A_i)
     true_sum+=ave;						//Vero valor medio					
     true_sum2+=ave*ave;
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
   int NoTot=1E4;
   int NoBlk=100;

   double S0=100.;
   double T=1.;
   double K=100.;
   double r=0.1;
   double vol=0.25;
   int N=100;
   double t=T/(double)(N);


   double meanval[NoTot];
   double St;
//----------------FIGURA 1--------------------------------

   char outfile1[20]="Es03.1.1.res";
   for(int i=0; i<NoTot; i++) 
   {
     St=S0*exp((r-0.5*pow(vol,2.))*T+vol*rnd.Gauss(0.,T));
     meanval[i]=exp(-r*T)*max(0.,St-K);
   }
   Blocks(NoTot, NoBlk, meanval, outfile1);

//----------------FIGURA 2--------------------------------

   char outfile2[20]="Es03.1.2.res";
   for(int i=0; i<NoTot; i++) 
   {
     St=S0;
     for(int ii=1; ii<=N; ii++)
       St=St*exp((r-0.5*pow(vol,2.))*(t*ii-t*(ii-1))+vol*rnd.Gauss(0.,1.)*sqrt(t*ii-t*(ii-1)));
     meanval[i]=exp(-r*T)*max(0.,St-K);
   }
   Blocks(NoTot, NoBlk, meanval, outfile2);

//----------------FIGURA 3--------------------------------

   char outfile3[20]="Es03.1.3.res";
   for(int i=0; i<NoTot; i++) 
   {
     St=S0*exp((r-0.5*pow(vol,2.))*T+vol*rnd.Gauss(0.,T));
     meanval[i]=exp(-r*T)*max(0.,K-St);
   }
   Blocks(NoTot, NoBlk, meanval, outfile3);

//----------------FIGURA 4--------------------------------

   char outfile4[20]="Es03.1.4.res";
   for(int i=0; i<NoTot; i++) 
   {
     St=S0;
     for(int ii=1; ii<=N; ii++)
       St=St*exp((r-0.5*pow(vol,2.))*(t*ii-t*(ii-1))+vol*rnd.Gauss(0.,1.)*sqrt(t*ii-t*(ii-1)));
     meanval[i]=exp(-r*T)*max(0.,K-St);
   }
   Blocks(NoTot, NoBlk, meanval, outfile4);


   rnd.SaveSeed();
   return 0;
}
