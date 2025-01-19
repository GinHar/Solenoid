#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include"matrix.h"
#include"integ.h"
#include"number_Gauss.h"

using namespace math;
using namespace std;



//Specific for a coil in the plane xy

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

  long double R=0.0079;//11cm
  long double I1=0.125;//10A
  long double mu=4*M_PI*pow(10,-7);
  long double L=0.704;//5m
  double N=3200;//

  double step=L/N;
  double L0=-L/2.0;
  
  //The positions of the coils are changed to calculate the effect of a solenoid

  matrix<long double> r2(3,1);
  matrix<long double> B(3,1);

  for(int i=0;i<N;i++){

    r2(0,0)=0;
    r2(1,0)=0;
    r2(2,0)=L0+i*step;

    B(0,0)=B(0,0)+Int_Gauss(0,1,R,r2,Bx_Integral);
    B(1,0)=B(1,0)+Int_Gauss(0,1,R,r2,By_Integral);
    B(2,0)=B(2,0)+Int_Gauss(0,1,R,r2,Bz_Integral);
    
  }

  B=mu*I1/4/M_PI*B;

  cout<<B<<endl;
  
}
