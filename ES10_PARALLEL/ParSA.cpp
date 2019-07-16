#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <cstring>
#include "ParSA.h"
#include "mpi.h"

int main(int argc, char* argv[]){

  MPI::Init(argc,argv);

  size = MPI::COMM_WORLD.Get_size();
  my_rank = MPI::COMM_WORLD.Get_rank();
  Inizialize(my_rank); 		//Generate the initial configuration
/*  if(!Check()){			//Checking if everithing is ok
    cerr << "Exiting with error -1 ..." << endl;
    return -1;
  }*/
  int cont;
  for(temp=temp_in, cont=0; cont<nstep; temp=temp/tau, cont++)
  {
    Mutation();
    Sort();					//Sorts the population from cost function
/*    if(!Check()){			//Checking if everithing is ok
      cerr << "Exiting with error -1 ..." << endl;
      return -1;
    }*/
    abs_min_val=CostFun(pop.front());
    MPI_Reduce(&abs_min_val, &abs_min_val_lec, 1, MPI::DOUBLE, MPI_MIN, 0, MPI::COMM_WORLD);
    PrintPartial(cont);		//Prints cost function of the best path and of mean value of the first half of the population
    if(my_rank==0) cerr << "Absolute best value: " << abs_min_val << endl;
  }

  PrintFinal();

  MPI::Finalize();
cerr << "TUTTO OK PER " << my_rank <<endl;
  return 0;
}


void Inizialize(int my_rank){//-----------------------------INITIALIZE--------------------------

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("/home/marco/LSN/PRNG/BASES/Primes");
   for(int i=0; i<my_rank; ++i) Primes >> p1 >> p2 ;
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

  for(int i=0; i<ncities; ++i)
    MPI_Bcast (coord[i].data(), coord[i].size(), MPI::DOUBLE, 0, MPI::COMM_WORLD);

  for(int i=0; i<ncities; ++i){
    for(int j=0; i<2; ++i)  cerr<<"rank: "<<my_rank<<", coord["<<i<<"]["<<j<<"]: "<<coord[i][j]<<endl;
  }
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

  ofstream Out, Position, Absbest, Lbest;
  char buffer1[30], buffer2[30], buffer3[30], buffer4[30];
  sprintf(buffer1, "conf.partial.%d", my_rank);
cerr << "My rank is " << my_rank << " and I'm printing on file " << buffer1 << endl;
  Out.open(buffer1);
  sprintf(buffer2, "position.%d", my_rank);
cerr << "My rank is " << my_rank << " and I'm printing on file " << buffer2 << endl;
  Position.open(buffer2);

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
  sprintf(buffer3, "output.lbest.%d", my_rank);
cerr << "My rank is " << my_rank << " and I'm printing on file " << buffer3 << endl;
  Lbest.open(buffer3, ios::app);
  Lbest << igeneration << setw(wd) << CostFun(pop.front()) << endl;
  Lbest.close();

  if(my_rank==0){
  //Printing absolute best value
    Absbest.open("output.absbest.0", ios::app);
    Absbest << igeneration << setw(wd) << abs_min_val_lec << endl;
    Absbest.close();
  }
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
