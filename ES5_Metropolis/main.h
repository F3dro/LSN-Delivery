#include <fstream>

//Global variables useful for Blocks method
double true_sum=0, true_sum2=0, sum=0;
int iCountBlk=0, BlockIndex=0;
double ave, sigma, true_ave, true_ave2;
std::ofstream outBlk;

int NoBlk=200;

int cont=0;
int const dim=3;
int nstep, state;
char distribution;
double x0[dim], x1[dim], passo, rapp;
double r=0;

int NoTot;

//Global Variables for generator
Random rnd;

class probabilita{
  public:
  probabilita(int dim){_dim=dim;};
  ~probabilita(){};

  double s1(double* x);
  double p2(double* x);
  private:
  double _r, _psi;
  int _dim;
};

int Input();
void Blocks(int , int , double , char* );
void InizializeRandom();
