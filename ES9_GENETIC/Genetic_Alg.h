#include <vector>
#include "/home/marco/LSN/PRNG/BASES/random.h"
using namespace std;

//Random numbers
int seed[4];
Random rnd;

vector<vector<int>> pop;	//Population vector of cities name vector
vector<vector<int>> new_pop;
vector<vector<double>> coord; 	//Coordinates in plan (2D) of the cities
vector<int> best;

int wd=15;



int mom, dad;
int ncities, nelements, ngenerations, conf_type, nhalf;
double _x1, _x2, _y1, _y2;
int n, m, st, el_dim; 		//Stuff for random mutations
double pm1, pm2, pm3, pm4, pm5, pcr;
double r;					//random probe number
int selected;

double R;				//Radius of the circumsphere on which the cities are placed

//functions
void Inizialize();
void Sort();
void Mutation();
bool Check();
void Crossover();
void NewGeneration();
void PrintPartial(int );
void PrintFinal();
void GeneratePosition(int );
double CostFun(vector<int>);
int Pbc(int);
int sign(double);
void Mutation2(int);
void Mutation3(int);
void Mutation4(int);
void Mutation5(int);
