#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include"res_sist.h"
#include"matrix.h"

using namespace math;
using namespace std;

//n es el numero de variables y ecuaciones, n/2 es el numero de puntos
//La variable x sus n primeros componentes son w y sus n últimos componentes son x

matrix<double> f(matrix<double> x){

  int n=x.RowNo();
  matrix<double> y(n,1);
  y.Null();
  
  for(int i=0;i<n;i++){

    for(int j=0;j<n/2.0;j++){

      if(i%2==0){
	y(i,0)=y(i,0)+x(j,0)*pow(x(j+n/2.0,0),i)-2/(i+1);
      } else {
	y(i,0)=y(i,0)+x(j,0)*pow(x(j+n/2.0,0),i);
      }
		
    }

  }

  return y;
  
}


matrix<double> fp(matrix<double> x){

  int n=x.RowNo();
  matrix<double> y(n,n);
  y.Null();

  for(int i=0;i<n;i++){

    for(int j=0;j<n;j++){
      
      if(j>=n/2.0){
	
	if(i==0){
	  y(i,j)=0;
	
	}else{
	  y(i,j)=i*x(j-n/2.0,0)*pow(x(j,0),i-1);  
	}
	
      }else{

	y(i,j)=pow(x(j,0),i);

      }
 
      
    }
  }

  return y;
  
}

void N(double n){


  
  matrix<double> x0(2*n,1);
  for(int i=0;i<2*n;i++){
    x0(i,0)=1;
  }
  cout<<x0<<endl;
  
  double eps=pow(10,-2);

  matrix<double> x=NR_sist(x0,f,fp,eps);

  cout<<x<<endl;
  
  
}
