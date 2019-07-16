
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>		//std::pow();
#include <iomanip>		//std::setw();
#include "../PRNG/BASES/random.h"
#include "main.h"

using namespace std;

int main (int argc, char *argv[]){

   int a=Input();
   if(a<0){
     cerr<<"Exiting with error "<<a<<endl;
     return -4;
   }
//Giving values to variables for blocks 
   NoTot=nstep;

   probabilita fzonda(dim);
   ofstream out;
   if(state==1)  out.open("output_GS.dat");
   else if(state==2) out.open("output_2p.dat");
//Start Metropolis with nstep points
   for(int j=0; j<nstep; j++){
     if(distribution=='U'){
       for(int i=0; i<dim; i++) x1[i]=x0[i]+rnd.Rannyu(-0.5,0.5)*passo; //Attemp step
     }
     else if(distribution=='G'){
       for(int i=0; i<dim; i++) x1[i]=x0[i]+rnd.Gauss(0,0.5)*passo; //Attemp step
     }
     if(state==2) rapp=fzonda.p2(x1)/fzonda.p2(x0);	//Using T=T^-1 so only p ratio matters
     else if(state==1) rapp=fzonda.s1(x1)/fzonda.s1(x0);	//Using T=T^-1 so only p ratio matters
     if(rnd.Rannyu() <= min(1., rapp) ){		//Accept-reject condition
       for(int i=0; i<dim; i++) x0[i]=x1[i];
       cont++;
     }
//Printing n+1 step on the output
     out<<j<<" ";
     for(int i=0; i<dim; i++){
       out<<x0[i]<<" ";
       r+=pow(x0[i],2.);
     }
     out<<endl;
     r=sqrt(r);
     if(state==2) Blocks(NoTot, NoBlk, r, "output_rmean_2p.dat");
     else if(state==1) Blocks(NoTot, NoBlk, r, "output_rmean_GS.dat");
     r=0;
   }
   out.close();
   cout<<"Rapporto di accettazione algoritmo Metropolis: "<< (double)cont/nstep<<endl;


   rnd.SaveSeed();
   return 0;
}




double probabilita::s1(double* x){
  _r=0;
  for(int i=0; i<_dim; i++) _r+=pow(x[i],2);
  _psi=exp(-sqrt(_r));
  return _psi*_psi;
}

double probabilita::p2(double* x){
  _r=0;
  for(int i=0; i<_dim; i++) _r+=pow(x[i],2);
//cerr<<"r^2 = "<<r_2<<endl;
  _psi=sqrt(_r)*exp(-sqrt(_r)/2.)*x[_dim-1]/sqrt(_r);
  return _psi*_psi;
}

int Input(){

   InizializeRandom();

//Reading input parameters
   ifstream in("input.dat");
   if(!in.good()){
     cerr<<"Error reading file: <input.dat> Exiting with error -1"<<endl;
     return -1;
   }
   in>>nstep;
   in>>passo;
   for(int i=0; i<dim; i++) in>>x0[i];
   in>>state;
   if(state!=1 && state!=2){
     cerr<<"Error in selected state: supported are 1s (state=1) or 2p (state=2)"<<endl;
     return -2;
   }
   in>>distribution;
   if(distribution!='U' && distribution!='G'){
     cerr<<"Error in selected distribution: supported are U (Uniform) and G (Gaussian)"<<endl;
     return -2;
   }
   in.close();

  return 0;
}

void Blocks(int NoTot, int NoBlk, double meanval, char* outfile)
{
//   NoTot=1E4;
//   NoBlk=100;
   int TrBlk=int(NoTot/NoBlk);
   iCountBlk++;
   if(iCountBlk<TrBlk){
     sum+=meanval;
//cerr<<"meanval = "<<meanval<<";	sum = "<<sum<<";	True_sum = "<<true_sum<<endl;
     return;
   }

   else if(iCountBlk==TrBlk){
     ave=(double)sum/(double)TrBlk;		//Media su ogni blocco (A_i)
     true_sum+=ave;						//Vero valor medio					
     true_sum2+=ave*ave;
     true_ave=true_sum/(double)(BlockIndex+1);
     true_ave2=true_sum2/(double)(BlockIndex+1);
     if((BlockIndex)==0){
       sigma=0;
       outBlk.open(outfile);		//First time rewrite the file
     }
     else{
       sigma=sqrt((double)(true_ave2-pow(true_ave,2))/BlockIndex);
       outBlk.open(outfile, ios::app);	//Other times append blocks results to the file
     }
     outBlk<<BlockIndex<<setw(15)<<true_ave<<setw(15)<<sigma<<endl;
 
     outBlk.close();
     iCountBlk=0;
     sum=0;
     BlockIndex++;
     return;
   }

   else{
   cerr<<"No matching condition, iCountBlk-value is "<<iCountBlk<<endl;
   return;
   }
}

void InizializeRandom(){
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
  return;
}
