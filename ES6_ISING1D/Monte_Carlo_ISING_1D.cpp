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
#include <ostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{

  system("rm output.mag.0 output.ene.0");
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> restart;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables


//initial configuration
  if(restart==1)
  {
    ifstream Readconf("config.final");
    for (int i=0; i<nspin; ++i) Readconf >> s[i];
    Readconf.close();
  }
  else
  {
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }

//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
      energy_old=Boltzmann(s[o], o);
      energy_new=Boltzmann(-1.*s[o], o);
      p=min(1.,exp(-(energy_new-energy_old)/temp));
      sm=rnd.Rannyu();
      if(sm<=p){
        s[o]*=-1;
        accepted++;
      }
      else s[o]*=1;
      attempted++;
    }
    else //Gibbs sampling
    {
      energy_up=Gibbs(+1, o, s);
      energy_down=Gibbs(-1, o, s);
      p=1./(1.+exp(+(energy_up-energy_down)/temp));  //probabilità che lo spin varrà +1
      sm=rnd.Rannyu();
      if(sm<=p) s[o]=+1;
      else s[o]=-1;
// INCLUDED MY CODE
    }
  }
}

double Boltzmann(int sm, int ip)  //sm valore dello spin prima di cambiare direzione (quindi prima temporalmente), ip intero che si e' mosso
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

double Gibbs(int sign, int pos, double *s){
  double sum1=0, sum2=0;
  s[pos]=sign;
  for(int i=0; i<nspin; ++i){
    sum1+=s[i]*s[Pbc(i+1)];
    sum2+=s[i]+s[Pbc(i+1)];
  }
  double ene = -J*sum1-0.5*h*sum2;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, c = 0.0, chi = 0.0, m=0.0, H;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     H = -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     u += H;
     chi += s[i] * s[i];
     m += s[i];
// INCLUDED MY CODE
  }
    for (int i=0; i<nspin; ++i)
  {
    c += H - u/(double)nspin;
  }
//cerr<<u<<endl;
  walker[iu] = u/(double)nspin;
  walker[ic] = beta*beta*c*c;
  walker[ix] = beta*chi;
  walker[im] = m;
// INCLUDEED MY CODE
}

void Reset(int iblk) //Reset block averages
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


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   ofstream EneT, HeatT, MagT, ChiT;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.0", ios::app);
    stima_u = blk_av[iu]/blk_norm; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    stima_c = blk_av[ic]/blk_norm; //Heat capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);

    stima_x = blk_av[ix]/blk_norm; //Susceptivity
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);

    Mag.open("output.mag.0", ios::app);
    stima_m = blk_av[im]/blk_norm; //magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Chi.close();

    if(iblk==nblk){
      EneT.open("output.ene.Temp",ios::app);
      HeatT.open("output.heat.Temp",ios::app);
      ChiT.open("output.chi.Temp",ios::app);
      MagT.open("output.mag.Temp",ios::app);
      EneT << setw(wd) << temp <<  setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
      HeatT << setw(wd) << temp <<  setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
      ChiT << setw(wd) << temp <<  setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
      MagT << setw(wd) << temp <<  setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
      EneT.close();
      HeatT.close();
      ChiT.close();
      MagT.close();
    }
// INCLUDE YOUR CODE HERE

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
