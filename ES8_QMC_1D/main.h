#include "../PRNG/BASES/random.h"

using namespace std;

//Global variables useful for Blocks method
double true_sum=0, true_sum2=0, sum=0;
int iCountBlk=0, BlockIndex=0;
double ave, sigma, true_ave, true_ave2;
ofstream outBlk;

//parameters, observables
const int m_props=100;
int n_props, iH;
double walker[m_props];
const int nbin=100;
const double xdx=3.0, xsx=-3.0, binsize=(xdx-xsx)/(double)nbin;
double histo[nbin];
const int wd=12;
int inhisto=0;

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima,err;

//variables
int nstep;
double x0, x1, passo, rapp, nblk;
double mu, sig;

void InizializeRandom();
int Input();
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move();
void Measure();
double potential(double x);
double psi(double x);
double psi_second(double x);
double Error(double,double,int);
void PrintFinal();

//Global Variables for generator
Random rnd;


