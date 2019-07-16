/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include <iomanip>		//setw();
#include "MolDyn_NVE.h"
#include "../PRNG/BASES/random.h"

using namespace std;

int main(){ 

  Input();             //Inizialization
  StartBlocks();
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure(istep);     //Properties measurement
        ConfXYZ(nconf);//Write actual configuration in XYZ format
        nconf += 1;
     }
  }
  ConfFinal();         //Write final configuration to restart
  DoBlocks();
  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;
  int restart;		

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

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
//-----------------NOW YOU CAN MODIFY------------------

  ReadInput.open("input.dat"); //Read input

  ReadInput >> restart;
  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> NoBlk;
  ReadInput >> material;

//For real simulations: it allows to perform quantitative simulations with real units
  if(material!=0){
    RealUnits=true;
    ifstream ReadMaterial;
    ReadMaterial.open("input_" + to_string(material) + ".material");
    ReadMaterial >> e_kb;
    ReadMaterial >> Sigma;
    ReadMaterial >> mass;
    ReadMaterial.close();
    Sigma*=1E-9;
  }

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  ip = 4; //Pressure
  n_props = 5; //Number of observables

  if(restart==0){
    nconf = 1;
//Read initial configuration
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();

//Prepare initial velocities
     cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
     double sumv[3] = {0.0, 0.0, 0.0};
     for (int i=0; i<npart; ++i){
       vx[i] = rnd.Rannyu() - 0.5;
       vy[i] = rnd.Rannyu() - 0.5;
       vz[i] = rnd.Rannyu() - 0.5;

       sumv[0] += vx[i];
       sumv[1] += vy[i];
       sumv[2] += vz[i];
     }
     for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
     double sumv2 = 0.0, fs;
     for (int i=0; i<npart; ++i){
       vx[i] = vx[i] - sumv[0];
       vy[i] = vy[i] - sumv[1];
       vz[i] = vz[i] - sumv[2];

       sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
     }
     sumv2 /= (double)npart;

     fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
cerr<<"fs = "<<fs<<endl;
     for (int i=0; i<npart; ++i){
       vx[i] *= fs;
       vy[i] *= fs;
       vz[i] *= fs;

       xold[i] = x[i] - vx[i] * delta;
       yold[i] = y[i] - vy[i] * delta;
       zold[i] = z[i] - vz[i] * delta;
     }
   }
   else if(restart==1){
//Read progressive number of configuration previously obtained from termalization operations
    ifstream ReadNconf;
    ReadNconf.open("nconf.txt");
    ReadNconf >> nconf;
    ReadNconf.close();

//Read initial configuration from previous simulation
    cout << "Read initial configuration from file config.final " << endl << endl;
    ReadConf.open("config.final");
    for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();

//Read initial old configuration from previous simulation
    cout << "Read initial old configuration from file config.old " << endl << endl;
    ReadConf.open("config.old");
    for (int i=0; i<npart; ++i){
      ReadConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
//cout<<x[i]<<"	"<<xold[i]<<endl;
    }
    ReadConf.close();

//Find the value of v(t)
    Move();
    double temp_att=0;

//Calcolo della temperatura attuale con metodo standard (usando velocità) //commentato
//    for(int i=0; i<npart; ++i){ temp_att+=(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i])/3.;}
//    temp_att/=npart;	

//Ricavo la temperatura attuale dalla media a blocchi della run precedente

    string temp_string_temp;
    if(RealUnits==true) temp_string_temp="output_ave_temp_" + to_string(material) + ".dat";
    else temp_string_temp="output_ave_temp.dat";
    ifstream ReadTemp(temp_string_temp);
    ReadTemp.clear();
    ReadTemp.seekg(-25, ReadTemp.end);
    ReadTemp>>temp_att;
    ReadTemp.close();

    if(RealUnits==true) temp_att /= e_kb;

    double fs = sqrt(temp / temp_att); 	// fs = trovo il fattore di scala dalla radice 
cerr<<"fs = "<<fs<<"	"<<temp_att<<"		"<<temp<<endl;
    for (int i=0; i<npart; ++i){	// del rapporto tra T_attesa e T_attuale
      vx[i] *= fs;			// e riscalo le velocita'
      vy[i] *= fs;
      vz[i] *= fs;
    }

    double appo;
    for (int i=0; i<npart; ++i){	//con le nuove velocita' trovo i passi iniziali e iniziali_old per ripartire
      appo= xold[i];
      xold[i] = Pbc( x[i] - 2.* delta * vx[i] );
      x[i]=appo;
      appo= yold[i];
      yold[i] = Pbc( y[i] - 2.* delta * vy[i] );
      y[i]=appo;
      appo= zold[i];
      zold[i] = Pbc( z[i] - 2.* delta * vz[i] );
      z[i]=appo;
    }
   }
   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
//cerr<<fx[i]<<"	"<<fy[i]<<"	"<<fz[i]<<endl;
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(int istep){ //Properties measurement
  int bin;
  double v, t, vij, p, pij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;

  if(RealUnits==true){
    Epot.open("output_epot_" + to_string(material) + ".dat",ios::app);
    Ekin.open("output_ekin_" + to_string(material) + ".dat",ios::app);
    Temp.open("output_temp_" + to_string(material) + ".dat",ios::app);
    Etot.open("output_etot_" + to_string(material) + ".dat",ios::app);
    Pres.open("output_pres_" + to_string(material) + ".dat",ios::app);
  }
  else{
    Epot.open("output_epot.dat",ios::app);
    Ekin.open("output_ekin.dat",ios::app);
    Temp.open("output_temp.dat",ios::app);
    Etot.open("output_etot.dat",ios::app);
    Pres.open("output_pres.dat",ios::app);
  }

  v = 0.0; //reset observables
  t = 0.0;
  p = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       pij = 48.0/pow(dr,12) - 24.0/pow(dr,6);

//Potential energy
       v += vij;
//Parzial pressure
       p += pij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total enery
    stima_pres = rho*stima_temp+(p/(double)npart)/(3.0*vol); //Pressure
//cerr<<stima_pot<<"		"<<stima_kin<<"		"<<stima_temp<<"		"<<stima_etot<<"		"<<stima_pres<<endl;
    if(RealUnits==true){
      stima_pot *= e_kb*kb;
      stima_kin *= e_kb*kb;
      stima_etot *= e_kb*kb;
      stima_pres *= e_kb*kb/pow(Sigma,3);
      stima_temp *= e_kb;
    }

//cerr<<stima_pot<<"		"<<stima_kin<<"		"<<stima_temp<<"		"<<stima_etot<<"		"<<stima_pres<<endl<<endl;

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Pres << stima_pres << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();


    ave_epot[istep/10 -1]=stima_pot;
    ave_ekin[istep/10 -1]=stima_kin;
    ave_temp[istep/10 -1]=stima_temp;
//cerr<<istep/10 -1<<") stima_temp: "<<stima_temp<<endl;
    ave_etot[istep/10 -1]=stima_etot;
    ave_pres[istep/10 -1]=stima_pres;

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

//Saving also the n-1 th cofiguration
  ofstream WriteConf_Old;

  cout << "Print final old configuration to file config.old " << endl << endl;
  WriteConf_Old.open("config.old");

  for (int i=0; i<npart; ++i){
    WriteConf_Old << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf_Old.close();

//Saving number of configuration
  ofstream WriteNconf;
  WriteNconf.open("nconf.txt");
  WriteNconf << nconf << endl;
  WriteNconf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;
  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
  return;
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

void StartBlocks(){//Creating vectors for block avarage
  NoTot=nstep/10.;
  ave_epot=new double[NoTot]; 
  ave_ekin=new double[NoTot]; 
  ave_etot=new double[NoTot]; 
  ave_temp=new double[NoTot]; 
  ave_pres=new double[NoTot]; 
  return;
}

void Blocks(double *meanval, string outfile)
{
   int TrBlk=int(NoTot/NoBlk);		//Throws per block
//cerr<<"TrBlk: "<<TrBlk<<"	NoBlk: "<<NoBlk<<"	NoTot: "<<NoTot<<"	file: "<<outfile<<endl;

   ifstream WriteNoBlk_count("NoBlk_count.txt");
   if(WriteNoBlk_count.good()) WriteNoBlk_count>>NoBlk_count;
   WriteNoBlk_count.close();

   ofstream out;
   out.open(outfile,ios::app);
   double ave;
   double sigma;
   double true_sum=0, true_sum2=0;
   double true_ave, true_ave2;
   int i;
   for(i=0; i<NoBlk; i++)
   {
     double sum=0;
     for(int j=0; j<TrBlk; j++)
        sum+=meanval[i*TrBlk+j];		//Integral to be computed
//cerr<<i*NoBlk+j<<") "<<meanval[i*TrBlk+j]<<endl;}
     ave=(double)sum/(double)TrBlk;		//Media su ogni blocco (A_i)
     true_sum+=ave;						//Vero valor medio					
     true_sum2+=ave*ave;
     true_ave=true_sum/(double)(i+1);
     true_ave2=true_sum2/(double)(i+1);
     if(i!=0)
       sigma=sqrt((double)(true_ave2-pow(true_ave,2))/i);
     else
       sigma=0;
     out<<i+NoBlk_count<<setw(15)<<true_ave<<setw(15)<<sigma<<endl;
   }
   NoBlk_count_appo=i;
   out.close();
   return;
}

void DoBlocks(){

  string output_ave_epot;
  string output_ave_etot;
  string output_ave_ekin;
  string output_ave_temp;
  string output_ave_pres;

  if(RealUnits==true){
    output_ave_epot="output_ave_epot_" + to_string(material) + ".dat";
    output_ave_etot="output_ave_etot_" + to_string(material) + ".dat";
    output_ave_ekin="output_ave_ekin_" + to_string(material) + ".dat";
    output_ave_temp="output_ave_temp_" + to_string(material) + ".dat";
    output_ave_pres="output_ave_pres_" + to_string(material) + ".dat";
  }

  else{
    output_ave_epot="output_ave_epot.dat";
    output_ave_etot="output_ave_etot.dat";
    output_ave_ekin="output_ave_ekin.dat";
    output_ave_temp="output_ave_temp.dat";
    output_ave_pres="output_ave_pres.dat";
  }

  cout<<"Printing average potenzial energy to file output_ave_epot.dat"<<endl;
  Blocks(ave_epot, output_ave_epot);

  cout<<"Printing average total energy to file output_ave_etot.dat"<<endl;
  Blocks(ave_etot, output_ave_etot);

  cout<<"Printing average kinetic energy to file output_ave_ekin.dat"<<endl;
  Blocks(ave_ekin, output_ave_ekin);

  cout<<"Printing average temperature to file output_ave_temp.dat"<<endl;
  Blocks(ave_temp, output_ave_temp);

  cout<<"Printing average pressure to file output_ave_pres.dat"<<endl;
  Blocks(ave_pres, output_ave_pres);

  ofstream WriteNoBlk_count("NoBlk_count.txt");
  WriteNoBlk_count << NoBlk_count+NoBlk_count_appo << endl;
  WriteNoBlk_count.close();

  delete[] ave_epot;
  delete[] ave_ekin;
  delete[] ave_etot;
  delete[] ave_temp;
  delete[] ave_pres;
  return;
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
