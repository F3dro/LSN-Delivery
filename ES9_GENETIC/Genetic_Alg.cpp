#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include "Genetic_Alg.h"

int main(){
  Inizialize(); 		//Generate the initial configuration
/*  if(!Check()){			//Checking if everithing is ok
    cerr << "Exiting with error -1 ..." << endl;
    return -1;
  }*/
  for(int i=0; i<ngenerations; ++i)
  {
    if(i%100==0) cout << "Generation " << i << " of " << ngenerations << endl;
    if(i>0 /*&& i%10==0*/) pop[pop.size()-1] = best;	//Replace last element with best element so far
    if(i==0) best = pop[0];
    Mutation();
    Sort();					//Sorts the population from cost function
/*    if(!Check()){			//Checking if everithing is ok
      cerr << "Exiting with error -1 ..." << endl;
      return -1;
    }*/
    PrintPartial(i);		//Prints cost function of the best path and of mean value of the first half of the population
    NewGeneration();	//Generates a new generation from the previous one
/*    if(!Check()){			//Checking if everithing is ok
      cerr << "Exiting with error -1 ..." << endl;
      return -1;
    }*/
  }

  PrintFinal();
  return 0;
}


void Inizialize(){//-----------------------------INITIALIZE--------------------------

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("/home/marco/LSN/PRNG/BASES/Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

// Read infos from file
  ifstream ReadInput;
  ReadInput.open("input.dat");
  ReadInput >> ncities;
  ReadInput >> ngenerations;
  ReadInput >> conf_type;
  ReadInput >> R;
// Reading probabilities for the mutations
  ReadInput >> pm1;		
  ReadInput >> pm2;
  ReadInput >> pm3;
  ReadInput >> pm4;
  ReadInput >> pm5;
  ReadInput >> pcr;
  ReadInput.close();

  nelements = ncities*ncities;

// Generate 2D positions of ncities cities
  GeneratePosition(conf_type);

// Generate initial population
  pop.resize(nelements);
  for(int i=0; i<nelements; ++i){
    for(int j=0; j<ncities; ++j){
      pop[i].push_back(j+1);
    }
    random_shuffle(pop[i].begin(), pop[i].end());
  }
  nhalf=pop.size()/2;
  return;
}

void Sort(){//--------------------------------SORT----------------------------------

  struct {
    bool operator()(vector<int> a, vector<int> b) const
    {   
        return CostFun(a) < CostFun(b);
    }   
  }shorter;

  sort(pop.begin(), pop.end(), shorter);

  if( CostFun(pop[0]) < CostFun(best) ) best = pop[0];

  return;
}

void PrintPartial(int igeneration){//-------------------PRINTPARTIAL----------------

  ofstream Out, Position, Lbest, Lave;
  Out.open("conf.partial");
  Position.open("position.0");

  for(int i=0; i<nelements; ++i){
    for(int j=0; j<ncities; ++j){
      Out << pop[i][j]; 		//Print name of cities ...
      if(i==0)
        //Here we save ordered position of cities for a sketch of the best path
        Position << coord[pop[i][j]-1][0] << setw(wd) << coord[pop[i][j]-1][1] << endl;
    }
    Out << '\t' << CostFun(pop[i]) << endl;		// ... and their cost function
  }
  Out.close();
  Position.close();

  //Printing best value, namely the first in the list, of Cost Function
  Lbest.open("output.lbest.0", ios::app);
  Lbest << igeneration << setw(wd) << CostFun(pop.front()) << endl;
  Lbest.close();

  //Printing average of Cost Function on the first half of the population
  Lave.open("output.lave.0", ios::app);
  double sum_best=0;
  for(int i=0; i<nhalf; ++i) sum_best+=CostFun(pop[i]);
  Lave << igeneration << setw(wd) << sum_best/nhalf << endl;
  Lave.close();
  return;
}

void Mutation(){//--------------------------------MUTATION--------------------------
  for(unsigned int i=0; i<pop.size()/2; ++i){
//Mutation are made more frequently on worse elements (with gaussian expectation)
    do{
      selected = int( -abs(rnd.Gauss(0., pop.size()/2)) + pop.size() );
    } while(selected < 0 || selected > pop.size());
//cerr << "Selected element for mutation: " << selected << " in a array of " << pop.size() << endl;
    r=rnd.Rannyu();
    if(r<pm1) Mutation2(selected);
    r=rnd.Rannyu();
    if(r<pm2) Mutation2(selected);
    r=rnd.Rannyu();
    if(r<pm3) Mutation3(selected);
    r=rnd.Rannyu();
    if(r<pm4) Mutation4(selected);
    r=rnd.Rannyu();
    if(r<pm5) Mutation5(selected);
  }
  return;

}

void Mutation1(int istep){
  auto first = pop[istep].begin();
  auto second = pop[istep].begin();
  advance(first,int(rnd.Rannyu(0,pop[istep].size())));
  advance(second,int(rnd.Rannyu(0,pop[istep].size())));
  iter_swap(first, second);
}

void Mutation2(int istep){
  n = int( rnd.Rannyu(0, pop[istep].size()) );
  rotate(pop[istep].rbegin(), pop[istep].rbegin() + n, pop[istep].rend());
  return;
}

void Mutation3(int istep){
  m = rnd.Rannyu(1, pop[istep].size()/2);
  n = rnd.Rannyu(1, pop[istep].size()/2);
  st = rnd.Rannyu(0, pop[istep].size()/2-n);
  for(int i=st; i<m; ++i){
    auto iter = pop[istep].begin();
    iter+=st;
    for(int j=0; j<n; ++j, ++iter){
      iter_swap(iter, iter+1);
    }
  }
  return;
}

void Mutation4(int istep){
  el_dim = pop[istep].size();
  n = int( rnd.Rannyu(1, el_dim/2) );
  m = int( rnd.Rannyu(0, el_dim-2*n) );
  st = int( rnd.Rannyu(0, el_dim-2*n-m) );
  swap_ranges(pop[istep].begin()+st, pop[istep].begin()+st+n, pop[istep].begin()+st+n+m);
  return;
}

void Mutation5(int istep){
  el_dim = pop[istep].size();
  n = int( rnd.Rannyu(1, el_dim) );
  st = int( rnd.Rannyu(0, el_dim-n) );
  reverse(pop[istep].begin()+st, pop[istep].begin()+st+n);
  return;
}

bool Check(){//--------------------------------CHECK------------------------------

  for(unsigned int i=0; i<pop.size(); ++i){
    for(unsigned int j=0; j<pop[i].size()-1; ++j){
      for(unsigned int jj=j+1; jj<pop[i].size(); ++jj){
        if(pop[i][j]==pop[i][jj]) return false;
      }
    }
  }

  return true;
}

void NewGeneration(){//----------------------------NEW_GENERATION---------------------

  for(int i=0; i<pop.size(); i+=2){
    do{
      mom = int( abs(rnd.Gauss(0., pop.size()/2)) );
    } while(mom < 0 || mom >= pop.size());
    do{
      dad = int( abs(rnd.Gauss(0., pop.size()/2)) );
    } while(dad < 0 || dad >= pop.size() || dad == mom);

    if(rnd.Rannyu()<pcr) Crossover();
    new_pop.push_back(pop[mom]);
    new_pop.push_back(pop[dad]);
  }
  pop=new_pop;
  new_pop.resize(0);
}

void Crossover(){//------------------------------CROSSOVER---------------------------
  n = int( rnd.Rannyu(1, ncities) );
  vector<int> mom_vec = pop[mom];
  vector<int> dad_vec = pop[dad];
  int cont_mom=0, cont_dad=0;
  bool pos_dad, pos_mom;

  for(int j=0; j<ncities /*|| (cont_dad==(ncities-n+1) && cont_mom==(ncities-n+1))*/; ++j){
    pos_dad=true;
    pos_mom=true;
    for(int i=0; i<n; ++i){
      if(dad_vec[j]==mom_vec[i]) pos_dad=false;
      if(mom_vec[j]==dad_vec[i]) pos_mom=false;
    }
    if(pos_dad){
      pop[mom][n+cont_mom]=dad_vec[j];
//cerr<<"dad sobstitution of "<< mom_vec[n+cont_mom] << " with " << dad_vec[j] << endl;
      cont_mom++;
    }
    if(pos_mom){
      pop[dad][n+cont_dad]=mom_vec[j];
//cerr<<"mom sobstitution of "<< dad_vec[n+cont_dad] << " with " << mom_vec[j] << endl;
      cont_dad++;
    }
  }

  return;
}

void GeneratePosition(int type){//---------------GENERATE_POSITION-------------------

coord.resize(ncities);
  if (type==1){
    for(int i=0; i<ncities; ++i){
      coord[i].push_back(rnd.Rannyu(-R,R));
      coord[i].push_back(sign(rnd.Rannyu(-1,1))*sqrt(R*R-coord[i][0]*coord[i][0]));
    }
  }

  else if (type==2) {
    for(int i=0; i<ncities; ++i){
      coord[i].push_back(rnd.Rannyu(-R,R));
      coord[i].push_back(rnd.Rannyu(-R,R));
    }
  }
  return;
}



void PrintFinal(){//-----------------------------PRINTFINAL--------------------------

  rnd.SaveSeed();

  return;
}

double CostFun(vector<int> element){//----------------COST_FUNCTION---------------
  double sum=0;
  for(unsigned int i=0; i<element.size(); ++i){
    _x1=coord[element[i]-1][0];
    _y1=coord[element[i]-1][1];
    _x2=coord[element[Pbc(i+1)]-1][0];
    _y2=coord[element[Pbc(i+1)]-1][1];
    sum+=fabs((_x2-_x1)*(_x2-_x1)+(_y2-_y1)*(_y2-_y1));
  }
  return sum;
}

//--------------------------------UTILITIES--------------------------------------

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= ncities) i = i - ncities;
    else if(i < 0) i = i + ncities;
    return i;
}

int sign(double x){
  if(x<0) return -1;
  else return +1;
}
