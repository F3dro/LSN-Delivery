
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>		//std::pow();
#include <iomanip>		//std::setw();
#include "../PRNG/BASES/random.h"

using namespace std;

double prob1(double* x, int dim){
  double r_2=0;
  for(int i=0; i<dim; i++) r_2+=pow(x[i],2);
  double psi=exp(-sqrt(r_2));
  return psi*psi;
}

double prob2(double x[], int dim){
  double r_2=0;
  for(int i=0; i<dim; i++) r_2+=pow(x[i],2);
cerr<<"r^2 = "<<r_2<<endl;
  double psi=sqrt(r_2)*exp(-sqrt(r_2)/2.)*x[dim-1]/sqrt(r_2);
/*
  double r_2=0;
  for(int i=0; i<dim; i++) r_2+=pow(x[i],2);
  return pow((1./8.)*sqrt(2./M_PI)sqrt(r_2)*exp(-sqrt(r_2)/2.)*x[2]/sqrt(r_2),2.);*/
  return psi*psi;
}

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
        sum+=meanval[i*TrBlk+j];		//Integral to be computed
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
   int NoBlk=200;

   int cont=0;
   int dim=3, nstep;
   double x0[dim], x1[dim], passo, rapp;
   double* r;

//Reading input parameters
   ifstream in("input.dat");
   if(!in.good()){
     cerr<<"Error reading file: <input.dat> Exiting with error -1"<<endl;
     return -1;
   }
   in>>nstep;
   in>>passo;
   for(int i=0; i<dim; i++) in>>x0[i];
   in.close();

//Defining variables for blocks 
   int NoTot=nstep;
   r=new double[nstep];
   for(int j=0; j<nstep; j++) r[j]=0;

   ofstream out("output_2p.dat");
//Start Metropolis with nstep points
   for(int j=0; j<nstep; j++){
     for(int i=0; i<dim; i++) x1[i]=x0[i]+rnd.Rannyu(-0.5,0.5)*passo; //Attemp step
     rapp=prob2(x1, dim)/prob2(x0, dim);		//Using T=T^-1 so only p ratio matters
     if(rnd.Rannyu() <= min(1., rapp) ){		//Accept-reject condition
       for(int i=0; i<dim; i++){
         x0[i]=x1[i];
       }
     cont++;
     }
//Printing n+1 step on the output
     out<<j<<" ";
     for(int i=0; i<dim; i++){
       out<<x0[i]<<" ";
       r[j]+=pow(x0[i],2.);
     }
     r[j]=sqrt(r[j]);
     out<<endl;
   }
   out.close();
   cout<<"Rapporto di accettazione algoritmo Metropolis: "<< (double)cont/nstep<<endl;

   Blocks(NoTot, NoBlk, r, "output_rmean_2p.dat");



   delete[] r;
   rnd.SaveSeed();
   return 0;
}

