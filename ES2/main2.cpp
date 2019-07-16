
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>		//std::pow();
#include <iomanip>		//std::setw();
#include "../PRNG/BASES/random.h"

using namespace std;
/*
class walk
{
  public:
    walk(): m_x(0.), m_y(0.), m_z(0.) {};
    walk(double x, double y, double z): m_x(x), m_y(y), m_z(z) {};
    void Addx(double val) {m_x+=val;}
    void Addy(double val) {m_y+=val;}
    void Addz(double val) {m_z+=val;}
    double Getx() {return m_x;}
    double Gety() {return m_y;}
    double Getz() {return m_z;}
    double GetMod2() {return pow(m_x,2.)+pow(m_y,2.)+pow(m_z,2.);}
    ~walk(){};
  private:
    double m_x, m_y, m_z;
} 
*/
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
   int nsteps=100;					//Steps of the random walk
//----------------PART 1----------------------------
   ofstream out("Es02.2.1.res");
   double ave[NoBlk][nsteps];
   double av2[NoBlk][nsteps];
   int sign=1;
   double passo=1.;
   double walk[3];
   double sum[nsteps];

   for(int i=0; i<NoBlk; i++)
   {
     for(int k=0; k<nsteps; k++)
     {
      sum[k]=0;
      ave[i][k]=0;
      av2[i][k]=0;
     }
     for(int j=0; j<TrBlk; j++)
      {
       walk[0]=0.;
       walk[1]=0.;
       walk[2]=0.;
       for(int k=0; k<nsteps;k++)
       {
        double x=rnd.Rannyu(-1.5,1.5);	//Generating random numbers between -1.5 and 1.5
//Get the direction of the step and move that way
        if(x>=0) sign = +1;
        else sign = -1;
        if(fabs(x)<0.5) walk[0]+=passo*sign;
        else if(fabs(x)<1.&&fabs(x)>0.5) walk[1]+=passo*sign;
        else if(fabs(x)<=1.5&&fabs(x)>1.) walk[2]+=passo*sign;       
        sum[k]+=pow(walk[0],2.)+pow(walk[1],2.)+pow(walk[2],2.);
       }
      }
      for(int k=0; k<nsteps; k++)
      {
       ave[i][k]=(double)sum[k]/(double)(TrBlk);
       av2[i][k]=pow(ave[i][k],2);
      }
   }

   double true_sum[NoBlk][nsteps];
   double true_sum2[NoBlk][nsteps];
   double sigma[NoBlk][nsteps];
  for(int k=0; k<nsteps; k++)
  {
   for(int i=0; i<NoBlk; i++)
   {
     true_sum[i][k]=0;
     true_sum2[i][k]=0;
     sigma[i][k]=0;
     for(int j=0; j<i+1; j++)
     {
       true_sum2[i][k]+=av2[j][k];
       true_sum[i][k]+=ave[j][k];
     }
     true_sum[i][k]/=(double)(i+1);
     true_sum2[i][k]/=(double)(i+1);
     if(i!=0)
       sigma[i][k]=sqrt((double)(true_sum2[i][k]-pow(true_sum[i][k],2))/((i+1)-1));
     else
       sigma[0][k]=0;
     if(i==NoBlk-1)
//     out<<k+1<<setw(15)<<sqrt(true_sum[i][k])<<setw(15)<<sqrt(sigma[i][k])<<endl;
     out<<k+1<<setw(15)<<sqrt(true_sum[i][k])<<setw(15)<<sigma[i][k]/(2.*sqrt(true_sum[i][k]))<<endl;

   }
  }
  
   out.close();



//----------------PART 2----------------------------
   ofstream out2("Es02.2.2.res");

   for(int i=0; i<NoBlk; i++)
   {
     for(int k=0; k<nsteps; k++)
     {
      sum[k]=0;
      ave[i][k]=0;
      av2[i][k]=0;
     }
     for(int j=0; j<TrBlk; j++)
      {
       walk[0]=0.;
       walk[1]=0.;
       walk[2]=0.;
       for(int k=0; k<nsteps;k++)
       {
        double phi=rnd.Rannyu(0.,2.*M_PI);		
        double theta=acos(1.-2.*rnd.Rannyu());	//Quantile of 1/2*sin(x)
        walk[0]+=passo*sin(theta)*cos(phi);
        walk[1]+=passo*sin(theta)*sin(phi);
        walk[2]+=passo*cos(theta);       
        sum[k]+=pow(walk[0],2.)+pow(walk[1],2.)+pow(walk[2],2.);
//if(i<1) cerr<<k<<setw(15)<<x<<setw(15)<<sum[k]<<endl;
       }
      }
      for(int k=0; k<nsteps; k++)
      {
       ave[i][k]=(double)sum[k]/(double)(TrBlk);
       av2[i][k]=pow(ave[i][k],2);
      }
   }


  for(int k=0; k<nsteps; k++)
  {
   for(int i=0; i<NoBlk; i++)
   {
     true_sum[i][k]=0;
     true_sum2[i][k]=0;
     sigma[i][k]=0;
     for(int j=0; j<i+1; j++)
     {
       true_sum2[i][k]+=av2[j][k];
       true_sum[i][k]+=ave[j][k];
     }
     true_sum[i][k]/=(double)(i+1);
     true_sum2[i][k]/=(double)(i+1);
     if(i!=0)
       sigma[i][k]=sqrt((double)(true_sum2[i][k]-pow(true_sum[i][k],2))/((i+1)-1));
     else
       sigma[0][k]=0;
//cerr<<k<<setw(15)<<i<<setw(15)<<true_sum[i][k]<<setw(15)<<true_sum2[i][k]<<setw(15)<<sigma[i][k]<<endl;
     if(i==NoBlk-1)
     out2<<k+1<<setw(15)<<sqrt(true_sum[i][k])<<setw(15)<<sigma[i][k]/(2.*sqrt(true_sum[i][k]))<<endl;

   }
  }
   out2.close();

   rnd.SaveSeed();
   return 0;
}
