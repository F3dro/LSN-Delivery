
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
   int NoBlk=100;
   int TrBlk=int(NoTot/NoBlk);		//Throws per block
//----------------PART 1----------------------------
   ofstream out("Es02.1.1.res");
   double ave[NoBlk];
   double av2[NoBlk];
   for(int i=0; i<NoBlk; i++)
   {
     double sum=0;
     for(int j=0; j<TrBlk; j++)
        sum+=M_PI_2*cos(M_PI_2*rnd.Rannyu());		//Integral to be computed
     ave[i]=(double)sum/(double)TrBlk;
     av2[i]=pow(ave[i],2);
   }
   double true_sum[NoBlk];
   double true_sum2[NoBlk];
   double sigma[NoBlk];
   for(int i=0; i<NoBlk; i++)
   {
     true_sum[i]=0;
     true_sum2[i]=0;
     for(int j=0; j<i+1; j++)
     {
       true_sum2[i]+=av2[j];
       true_sum[i]+=ave[j];
     }
     true_sum[i]/=(double)(i+1);
     true_sum2[i]/=(double)(i+1);
     if(i!=0)
       sigma[i]=sqrt((double)(true_sum2[i]-pow(true_sum[i],2))/((i+1)-1));
     else
       sigma[0]=0;
     out<<i<<setw(15)<<true_sum[i]<<setw(15)<<sigma[i]<<endl;
   }
   out.close();

//----------------PART 2----------------------------
   ofstream out2("Es02.1.2.res");
   for(int i=0; i<NoBlk; i++)
   {
     double sum=0;
     ave[i]=0;
     av2[i]=0;
     for(int j=0; j<TrBlk; j++)
      {
      double x=+1.-sqrt(1-rnd.Rannyu());		//Generating points distributed with f(x)=2x-x^2
      sum+=M_PI_2*cos(M_PI_2*x)/(2.*(1.-x));
      }
     ave[i]=(double)sum/(double)(TrBlk);
     av2[i]=pow(ave[i],2);
   }

   for(int i=0; i<NoBlk; i++)
   {
     true_sum[i]=0;
     true_sum2[i]=0;
     sigma[i]=0;
     for(int j=0; j<i+1; j++)
     {
       true_sum2[i]+=av2[j];
       true_sum[i]+=ave[j];
     }
     true_sum[i]/=(double)(i+1);
     true_sum2[i]/=(double)(i+1);
//cerr<<true_sum[i]<<"	"<<true_sum2[i]<<endl;
     if(i!=0)
       sigma[i]=sqrt((double)(true_sum2[i]-pow(true_sum[i],2))/((i+1)-1));
     else
       sigma[0]=0;
     out2<<i<<setw(15)<<true_sum[i]<<setw(15)<<sigma[i]<<endl;
   }
   out2.close();


   rnd.SaveSeed();
   return 0;
}
