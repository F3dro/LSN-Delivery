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
   int NoBlk=100;
   int TrBlk=int(NoTot/NoBlk);
//----------------PART 1----------------------------
   ofstream out("Es01.1.1.res");
   double ave[NoBlk];
   double av2[NoBlk];
   for(int i=0; i<NoBlk; i++)
   {
     double sum=0;
     for(int j=0; j<TrBlk; j++)
        sum+=rnd.Rannyu();
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
//cerr<<true_sum[i]<<"	"<<true_sum2[i]<<endl;
     if(i!=0)
       sigma[i]=sqrt((double)(true_sum2[i]-pow(true_sum[i],2))/((i+1)-1));
     else
       sigma[0]=0;
     out<<i<<setw(15)<<true_sum[i]<<setw(15)<<sigma[i]<<endl;
   }
   out.close();

//----------------PART 2----------------------------
   ofstream out2("Es01.1.2.res");
   for(int i=0; i<NoBlk; i++)
   {
     double sum=0;
     ave[i]=0;
     av2[i]=0;
     for(int j=0; j<TrBlk; j++)
        sum+=pow((rnd.Rannyu()-0.5),2);
     ave[i]=(double)sum/(double)TrBlk;
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

/*
//-----------------TRYING----------------------
   Primes.open("../PRNG/BASES/Primes");
   Primes.clear();
   Primes.close();
   Primes.open("../PRNG/BASES/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   input.open("../PRNG/BASES/seed.in");
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


   ofstream out3("Es01.1.2_2.res");
   for(int i=0; i<NoBlk; i++)
   {
     double sum=0;
     ave[i]=0;
     av2[i]=0;
     for(int j=0; j<TrBlk; j++)
        sum+=pow((rnd.Rannyu()-0.5),2);
     ave[i]=(double)sum/(double)TrBlk;
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
     out3<<true_sum[i]<<setw(15)<<sigma[i]<<endl;
   }
   out3.close();
*/
//----------------------PART 3--------------------------
   int M=100;		//Number of intervals
   int n=1E4;		//Number of throws
   double chi2=0, ni[n];
   int cont=0;
   ofstream out3("Es01.1.3.res");
   for(int k=0; k<100; k++)
   {
     chi2=0;
     for(int j=0; j<n; j++) ni[j]=rnd.Rannyu();
     for(int i=0; i<M; i++)
     {
       cont=0;
       for(int j=0; j<n; j++) 
         {if(ni[j]>=((double)i/(double)(M))&&ni[j]<((double)(i+1)/(double)(M))) cont++;}
       chi2+=(double)pow((cont-n/M),2.)/(double)(n/M);
     }
     out3<<k<<setw(15)<<chi2<<endl;
   }
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
