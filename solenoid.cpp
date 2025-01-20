#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include"matrix.h"
#include"integ.h"

using namespace math;
using namespace std;



//Specific for a coil in the plane XY
//Functions that need to be integrated to calculated the effect of one coil
//Each function calculates the magnetic field in one axis.

double Bx_Integral(long double R,long double t, matrix<long double> r2){

  matrix<long double> r(3,1);
  r(0,0)=r2(0,0)-R*cos(2*M_PI*t);
  r(1,0)=r2(1,0)-R*sin(2*M_PI*t);
  r(2,0)=r2(2,0);
  long double B=R*2*M_PI*cos(2*M_PI*t)*r(2,0)/pow(r.Norm(),3);
  return B;
  
}

double By_Integral(long double R,long double t, matrix<long double> r2){

  matrix<long double> r(3,1);
  r(0,0)=r2(0,0)-R*cos(2*M_PI*t);
  r(1,0)=r2(1,0)-R*sin(2*M_PI*t);
  r(2,0)=r2(2,0);
  long double B=R*2*M_PI*sin(2*M_PI*t)*r(0,0)/pow(r.Norm(),3);
  return B;
  
}

double Bz_Integral(long double R,long double t, matrix<long double> r2){

  matrix<long double> r(3,1);
  r(0,0)=r2(0,0)-R*cos(2*M_PI*t);
  r(1,0)=r2(1,0)-R*sin(2*M_PI*t);
  r(2,0)=r2(2,0);
  long double B=-R*2*M_PI*sin(2*M_PI*t)*r(1,0)/pow(r.Norm(),3)-R*2*M_PI*cos(2*M_PI*t)*r(0,0)/pow(r.Norm(),3);
  return B;
  
}

int main(){

  //Data
  long double R=0.005;//0.0079m
  long double I1=1;//0.125A
  long double mu=4*M_PI*pow(10,-7);
  long double L=0.5;//0.704m
  double N=320;//

  double step=L/N;
  
  //The positions of the coils are changed to calculate the effect of a solenoid

  matrix<long double> r2(3,1);
  matrix<long double> B(3,1); //It is a variable that is going to be changed
  
  
  // Data
  ofstream outfile("Field_XZ.txt");


  //We calculated the magnetic field in the XZ plane

  double l=1; // 1m² is the region that we are going to see
  double space=0.02; // The distance between each dot is 0.1m
  double n=l/space; // The number of dots in each axis
  
  for(int i=0;i<n;i++){
    
    for(int j=0;j<n;j++){
      
      for(int k=0;k<N;k++){

	r2(0,0)=i*space-l/2.0;
	r2(2,0)=k*step-L/2.0+j*space-l/2.0;

	B(0,0)=B(0,0)+Int_Gauss(0,1,R,r2,Bx_Integral);
	B(2,0)=B(2,0)+Int_Gauss(0,1,R,r2,Bz_Integral);

      }
      B=mu*I1/4/M_PI*B;
      outfile <<r2(0,0)<<" "<<j*space-l/2.0<<" "<<B(0,0)<<" "<<B(2,0)<<endl;
      B(0,0)=0;
      B(2,0)=0;
              
    }
  }

  outfile.close();

  cout<<"Magnetic field calculated"<<endl;
  
}

