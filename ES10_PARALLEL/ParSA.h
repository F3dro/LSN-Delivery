#include <vector>
#include "/home/marco/LSN/PRNG/BASES/random.h"
using namespace std;

//Random numbers
int seed[4];
Random rnd;

//Parallel stuff
int my_rank, size;
double abs_min_val, abs_min_val_lec;

vector<vector<int>> pop;	//Population vector of cities name vector
vector<vector<int>> new_pop;
vector<vector<double>> coord; 	//Coordinates in plan (2D) of the cities
vector<int> appo;

int wd=15;

//Simulated Annealing
double temp, temp_in, temp_fin, T_step;
double tau;
int nstep_per_T, nstep;
int acc1, acc3, acc4, acc5, att1, att3, att4, att5;

int ncities, nelements, conf_type, nhalf;
double _x1, _x2, _y1, _y2;
int n, m, st, el_dim; 		//Stuff for random mutations
double r, p;					//random probe number
int selected;

double R;				//Radius of the circumsphere on which the cities are placed

//functions
void Inizialize(int);
void Sort();
void Mutation();
bool Check();
void PrintPartial(int );
void PrintFinal();
void GeneratePosition(int );
double CostFun(vector<int>);
int Pbc(int);
int sign(double);
void Mutation1(int);
void Mutation3(int);
void Mutation4(int);
void Mutation5(int);
