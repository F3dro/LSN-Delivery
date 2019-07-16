#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include "SA.h"

int main(){
  Inizialize(); 		//Generate the initial configuration
/*  if(!Check()){			//Checking if everithing is ok
    cerr << "Exiting with error -1 ..." << endl;
    return -1;
  }*/
  int cont;
  for(temp=temp_in, cont=0; cont<nstep; temp=temp/tau, cont++)
  {
    if(cont%100==0) cout << "Step " << cont << " of " << nstep << endl;
    Mutation();
    Sort();					//Sorts the population from cost function
/*    if(!Check()){			//Checking if everithing is ok
      cerr << "Exiting with error -1 ..." << endl;
      return -1;
    }*/
    PrintPartial(cont);		//Prints cost function of the best path and of mean value of the first half of the population
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
  ReadInput >> nstep;
  ReadInput >> conf_type;
  ReadInput >> R;
// Reading probabilities for the mutations
  ReadInput >> temp_in;
  ReadInput >> tau;
  ReadInput >> nstep_per_T;
  ReadInput.close();

  nelements = 1;

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

  return;
}

void PrintPartial(int igeneration){//-------------------PRINTPARTIAL----------------

  ofstream Out, Position, Lbest, Temp;
  Out.open("conf.partial.one");
  Position.open("position.one");

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
  Lbest.open("output.lbest.one", ios::app);
  Lbest << igeneration << setw(wd) << CostFun(pop.front()) << endl;
  Lbest.close();

  //Printing temperature trend
  Temp.open("output.temp.one", ios::app);
  Temp << igeneration << setw(wd) << temp << endl;
  Temp.close();
  return;
}

void Mutation(){//--------------------------------MUTATION--------------------------
  att1=0, att3=0, att4=0, att5=0;
  acc1=0, acc3=0, acc4=0, acc5=0;
  for(int i=0; i<nstep_per_T; ++i){
    selected=0;
//cerr << "Selected element for mutation: " << selected << " in a array of " << pop.size() << endl;
    Mutation1(selected);
    Mutation3(selected);
    Mutation4(selected);
    Mutation5(selected);
  }
  cout << "T = " << temp << ", accettazione per la mutazione 1 : " << (double)acc1/att1 << endl;
  cout << "T = " << temp << ", accettazione per la mutazione 3: " << (double)acc3/att3 << endl;
  cout << "T = " << temp << ", accettazione per la mutazione 4: " << (double)acc4/att4 << endl;
  cout << "T = " << temp << ", accettazione per la mutazione 5: " << (double)acc5/att5 << endl;
cerr << endl;
  return;

}

void Mutation1(int istep){
  appo=pop[istep];
  auto first = pop[istep].begin();
  auto second = pop[istep].begin();
  advance(first,int(rnd.Rannyu(0,pop[istep].size())));
  advance(second,int(rnd.Rannyu(0,pop[istep].size())));
  iter_swap(first, second);
  p=min(1., exp(-( CostFun(pop[istep]) - CostFun(appo) )/temp));
  r=rnd.Rannyu();
  if(r<=p) acc1++;
  else pop[istep]=appo;
  att1++;
}

void Mutation3(int istep){
  appo=pop[istep];
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
  p=min(1., exp(-( CostFun(pop[istep]) - CostFun(appo) )/temp));
  r=rnd.Rannyu();
  if(r<=p) acc3++;
  else pop[istep]=appo;
  att3++;
  return;
}

void Mutation4(int istep){
  appo=pop[istep];
  el_dim = pop[istep].size();
  n = int( rnd.Rannyu(1, el_dim/2) );
  m = int( rnd.Rannyu(0, el_dim-2*n) );
  st = int( rnd.Rannyu(0, el_dim-2*n-m) );
  swap_ranges(pop[istep].begin()+st, pop[istep].begin()+st+n, pop[istep].begin()+st+n+m);
  p=min(1., exp(-( CostFun(pop[istep]) - CostFun(appo) )/temp));
  r=rnd.Rannyu();
  if(r<=p) acc4++;
  else pop[istep]=appo;
  att4++;
  return;
}

void Mutation5(int istep){
  appo=pop[istep];
  el_dim = pop[istep].size();
  n = int( rnd.Rannyu(1, el_dim) );
  st = int( rnd.Rannyu(0, el_dim-n) );
  reverse(pop[istep].begin()+st, pop[istep].begin()+st+n);
  p=min(1., exp(-( CostFun(pop[istep]) - CostFun(appo) )/temp));
  r=rnd.Rannyu();
  if(r<=p) acc5++;
  else pop[istep]=appo;
  att5++;
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
