
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>		//std::pow();
#include <iomanip>		//std::setw();
#include "../PRNG/BASES/random.h"

using namespace std;

void Blocks(int NoTot, int NoBlk, double *meanval,char* outfile)
{
//   int NoTot=1E4;
//   int NoBlk=100;
   int TrBlk=int(NoTot/NoBlk);		//Throws per block
//----------------PART 1----------------------------
   ofstream out(outfile);
   double ave[NoBlk];
   double av2[NoBlk];
   for(int i=0; i<NoBlk; i++)
   {
     double sum=0;
     for(int j=0; j<TrBlk; j++)
        sum+=meanval[i*NoBlk+j];		//Integral to be computed
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
   int TrBlk=int(NoTot/NoBlk);		//Throws per block

   double S0=100.;
   double T=1.;
   double K=100.;
   double r=0.1;
   double vol=0.25;
   int N=100;
   double t=T/(double)(N);
//----------------PART 1----------------------------
   ofstream out("Es03.1.1_backup.res");
   double ave[NoBlk];
   double av2[NoBlk];
   double ave_2[NoBlk];
   double av2_2[NoBlk];
   for(int i=0; i<NoBlk; i++)
   {
     double sum=0;
     double sum_2=0;
     for(int j=0; j<TrBlk; j++)
     {
        double St=S0*exp((r-0.5*pow(vol,2.))*T+vol*rnd.Gauss(0.,T));
        sum+=exp(-r*T)*max(0.,St-K);		//Integral to be computed
        St=S0;
        for(int ii=1; ii<=N; ii++)
         St=St*exp((r-0.5*pow(vol,2.))*(t*ii-t*(ii-1))+vol*rnd.Gauss(0.,1.)*sqrt(t*ii-t*(ii-1)));
        sum_2+=exp(-r*T)*max(0.,St-K);
     }
     ave[i]=(double)sum/(double)TrBlk;
     av2[i]=pow(ave[i],2);
     ave_2[i]=(double)sum_2/(double)TrBlk;
     av2_2[i]=pow(ave_2[i],2);
   }
   double true_sum[NoBlk];
   double true_sum2[NoBlk];
   double sigma[NoBlk];
   double true_sum_2[NoBlk];
   double true_sum2_2[NoBlk];
   double sigma_2[NoBlk];
   for(int i=0; i<NoBlk; i++)
   {
     true_sum[i]=0;
     true_sum2[i]=0;
     true_sum_2[i]=0;
     true_sum2_2[i]=0;
     for(int j=0; j<i+1; j++)
     {
       true_sum2[i]+=av2[j];
       true_sum[i]+=ave[j];
       true_sum2_2[i]+=av2_2[j];
       true_sum_2[i]+=ave_2[j];
     }
     true_sum[i]/=(double)(i+1);
     true_sum2[i]/=(double)(i+1);
     true_sum_2[i]/=(double)(i+1);
     true_sum2_2[i]/=(double)(i+1);
     if(i!=0)
     {
       sigma[i]=sqrt((double)(true_sum2[i]-pow(true_sum[i],2))/((i+1)-1));
       sigma_2[i]=sqrt((double)(true_sum2_2[i]-pow(true_sum_2[i],2))/((i+1)-1));
     }
     else
     {
       sigma[0]=0;
       sigma_2[0]=0;
     }
     out<<i<<setw(15)<<true_sum[i]<<setw(15)<<sigma[i]<<setw(25)<<true_sum_2[i]<<setw(15)<<sigma_2[i]<<endl;
   }
   out.close();

//----------------PART 2----------------------------
   ofstream out2("Es03.1.2_backup.res");

   for(int i=0; i<NoBlk; i++)
   {
     double sum=0;
     double sum_2=0;
     for(int j=0; j<TrBlk; j++)
     {
        double St=S0*exp((r-0.5*pow(vol,2.))*T+vol*rnd.Gauss(0.,T));
        sum+=exp(-r*T)*max(0.,K-St);		//Integral to be computed
        St=S0;
        for(int ii=1; ii<=N; ii++)
         St=St*exp((r-0.5*pow(vol,2.))*(t*ii-t*(ii-1))+vol*rnd.Gauss(0.,1.)*sqrt(t*ii-t*(ii-1)));
        sum_2+=exp(-r*T)*max(0.,K-St);
     }
     ave[i]=(double)sum/(double)TrBlk;
     av2[i]=pow(ave[i],2);
     ave_2[i]=(double)sum_2/(double)TrBlk;
     av2_2[i]=pow(ave_2[i],2);
   }

   for(int i=0; i<NoBlk; i++)
   {
     true_sum[i]=0;
     true_sum2[i]=0;
     true_sum_2[i]=0;
     true_sum2_2[i]=0;
     for(int j=0; j<i+1; j++)
     {
       true_sum2[i]+=av2[j];
       true_sum[i]+=ave[j];
       true_sum2_2[i]+=av2_2[j];
       true_sum_2[i]+=ave_2[j];
     }
     true_sum[i]/=(double)(i+1);
     true_sum2[i]/=(double)(i+1);
     true_sum_2[i]/=(double)(i+1);
     true_sum2_2[i]/=(double)(i+1);
     if(i!=0)
     {
       sigma[i]=sqrt((double)(true_sum2[i]-pow(true_sum[i],2))/((i+1)-1));
       sigma_2[i]=sqrt((double)(true_sum2_2[i]-pow(true_sum_2[i],2))/((i+1)-1));
     }
     else
     {
       sigma[0]=0;
       sigma_2[0]=0;
     }
     out2<<i<<setw(15)<<true_sum[i]<<setw(15)<<sigma[i]<<setw(25)<<true_sum_2[i]<<setw(15)<<sigma_2[i]<<endl;
   }
   out2.close();


   rnd.SaveSeed();
   return 0;
}
