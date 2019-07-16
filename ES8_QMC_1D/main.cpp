#include <iostream>
#include <fstream>
#include <string>
#include <cmath>		//std::pow();
#include <iomanip>		//std::setw();
#include "main.h"








int main (int argc, char *argv[]){

   InizializeRandom();
   if(Input()==-1) return -1;

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
//Start Metropolis with nstep points per block
    for(int j=0; j<nstep; j++){
     Move();
     Measure();
     Accumulate();
    }
  Averages(iblk);
  } 

  PrintFinal();

  rnd.SaveSeed();
  return 0;
}







int Input(){//------------------------------------INPUT----------------------------
//Reading input parameters
  ifstream in("input.dat");
  if(!in.good()){
    cerr<<"Error reading file: <input.dat> Exiting with error -1"<<endl;
    return -1;
  }
  in >> nstep >> passo >> nblk >> mu >> sig;
  in>>x0;
  in.close();

//Prepare arrays for measurements
  iH = 0; //energy
  n_props = 1; //Number of observables

  for(int i=0; i<nbin; ++i) histo[i]=0.;

  return 0;
}

void Move(){//------------------------------------MOVE----------------------------
  x1=x0+rnd.Rannyu(-0.5,0.5)*passo; //Attemp step
  rapp=(psi(x1)*psi(x1))/(psi(x0)*psi(x0));	//Using T=T^-1 so only p ratio matters
  if(rnd.Rannyu() <= min(1., rapp) ){		//Accept-reject condition
    x0=x1;
    accepted++;
    inhisto++;
    int pos=floor(x1/binsize)+nbin/2;
    histo[pos]++;
  }
  attempted++;
  return;
}

void Measure(){
  walker[iH]=(-0.5*psi_second(x0)+potential(x0)*psi(x0))/psi(x0);
  return;
}

double potential(double x){
  return pow(x,4)-5.0*x*x/2.0;
}

double psi(double x){
  return exp( -((x-mu)*(x-mu))/(2.0*sig*sig) ) + exp( -((x+mu)*(x+mu))/(2.0*sig*sig) );
}

double psi_second(double x){
  return exp( -((x+mu)*(x+mu))/(2.0*sig*sig) ) * ( -sig*sig + exp( (2.0*mu*x)/(sig*sig) )*(-mu-sig+x)*(-mu+sig+x) +(mu+x)*(mu+x)  ) / pow(sig,4.0);
//  return - exp(-(x-mu)*(x-mu)/(2.0*sig*sig))/sig/sig - exp(-(x+mu)*(x+mu)/(2.0*sig*sig))/sig/sig + (x-mu)*(x-mu)*exp(-(x-mu)*(x-mu)/(2.0*sig*sig))/pow(sig,4) + (x+mu)*(x+mu)*exp(-(x+mu)*(x+mu)/(2.0*sig*sig))/pow(sig,4);
}

void Reset(int iblk) //--------------------Reset block averages-----------------------
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Accumulate(void) //------------------Update block averages------------------------
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //----------------Print results for current block-------------
{
   ofstream out;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    out.open("output.ene.0",ios::app);
    
    stima = blk_av[iH]/blk_norm; //(H*psi)/psi
    glob_av[iH] += stima;
    glob_av2[iH] += stima*stima;
    err=Error(glob_av[iH],glob_av2[iH],iblk);

//Potential energy per particle
    out << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[iH]/(double)iblk << setw(wd) << err << setw(wd) << mu << setw(wd) << sig << endl;


    cout << "----------------------------" << endl << endl;

    out.close();
}

void PrintFinal(){//----------------------------PRINT FINAL---------------------------
  ofstream out_histo;
  out_histo.open("output.psi2.0");
  for(int i=0; i<nbin; ++i) 
    out_histo << xsx+i*binsize << setw(wd) << (double)histo[i]/(double)inhisto << endl;
  out_histo.close();
  return;
}

void InizializeRandom(){//-------------------------RANDOM STUFF-----------------------
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

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}
 


