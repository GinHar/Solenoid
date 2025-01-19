#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include"matrix.h"

using namespace math;
using namespace std;

double Int_Trap(double x0, double xf, int n, double &Err, double(*y)(double), double(*yp)(double)){

  double tau=(xf-x0)/n;
  matrix<double> x(n+1,1);
  x(0,0)=x0;
  
  for(int i=1;i<n;i++){

    x(i,0)=x(i-1,0)+tau;

  }

  x(n,0)=xf;

  double val=0;
  
  for(int i=1;i<n;i++){

    val=val+y(x(i,0));

  }

  val=0.5*y(x(0,0))+val+0.5*y(x(n,0));
  val=tau*val;

  Err=(tau*tau/12.0)*(yp(xf)-yp(x0));

  return val;

}

double Int_Simpson(double x0, double xf, int n, double &Err, double(*y)(double), double(*yppp)(double)){

  if(n%2==0){
    
    double val=0;
    double tau=(xf-x0)/n;
    matrix<double> x(n+1,1);
    x(0,0)=x0;
  
    for(int i=1;i<n;i++){

      x(i,0)=x(i-1,0)+tau;

    }

    x(n,0)=xf;


    for(int i=1;i<n;i++){

      if(i%2==0){

	val=val+2*y(x(i,0));
	
      }else{

	val=val+4*y(x(i,0));
	
      }

    }

    val=y(x(0,0))+val+y(x(n,0));
    val=tau/3.0*val;

    Err=-1.0/180.0*pow(tau,4)*(yppp(xf)-yppp(x0));

    return val;
    

  }else{

    cout<<"El numero de intervalos no son pares"<<endl;
    double val=0;
    return val;
    
  }

}




long double Int_Gauss(long double a,long double b,long double R, matrix<long double> r2,double(*y)(long double,long double, matrix<long double>)){

  int n;
  long double w;
  long double x;
  
  ifstream in("Gauss_100n.txt");

  long double val=0;

  for(int i=0;i<100;i++){

    in >> n;
    in >> x;
    in >> w;
    val=val+w*y(R,((b+a)-(b-a)*x)/2,r2);

  }

  return (b-a)/2*val;
  
}
