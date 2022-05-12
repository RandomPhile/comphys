#include <iostream>
#include <armadillo>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace arma;
typedef Mat<double> mat;
typedef Col<double> vec;

const int n=5; //particelle a lato del cubo
const int N = n*n*n;
const double pi = 3.141592653589793;
const double rho = 0.01;
mat r = zeros(N,3);
mat vR = zeros(N-1,3);
vec R = zeros(N-1);
vec R6 = zeros(N-1);
vec R12 = zeros(N-1);
const double T = 1.1;
const double sigma_v = sqrt(T);
const double dl = pow(1./rho,1./3.);
const double L = n*dl;
const double delta = dl;

int MRT(mat r, double delta);
double pot(mat r);

double cont = 0;

int main(void)
{
	int step = 1000;
	vec U = zeros(step);
	
	
	//inizializziamo posizioni
	for (int k = 0; k<n; ++k)
    	{ 
        	for(int l = 0; l<n; ++l)
        	{
		    	for(int m = 0; m<n; ++m)
		    	{
		        int ii = (k*n*n+l*n+m);
		        r(ii,0) = k*dl;
		        r(ii,1) = l*dl;
		        r(ii,2) = m*dl;
		    	} 
        	}
	}
	
	U(0) = pot(r);
	for (int kk = 1; kk <step; ++kk)
	{
		MRT(r,delta);
		U(kk) = pot(r);
	}

	U.save("U.txt",csv_ascii);
	r.save("rx.txt",csv_ascii);
	
	cout << cont <<endl;
	return 0;
}


double pot(mat r)
{
	vec UU = zeros(N);
	double U;
	for (int ii = 0; ii < N; ++ii)
	{
		for (int jj = 0; jj<N && jj != ii; ++jj)
		{
			vR.row(jj) = r.row(ii) - r.row(jj)  ;
			vR.row(jj) = vR.row(jj) - L*round(vR.row(jj)/L) ;
			
			R(jj) = sqrt(  accu( vR.row(jj)%vR.row(jj) )  );	
		}
		
		R6 = square(R)%square(R)%square(R);
		R12 = R6%R6;
		UU(ii) = accu(4/R12-4/R6));
			
	} 
	U = accu(UU);
	return U;	
}


int MRT(mat r, double delta)
{
	double U,U_,A;
	vec dr = zeros(3);
	
	for (int ii = 0; ii < N; ++ii)
	{
		
		for (int jj = 0; jj<N && jj != ii; ++jj)
		{
			vR.row(jj) = r.row(ii) - r.row(jj)  ;
			vR.row(jj) = vR.row(jj) - L*round(vR.row(jj)/L) ;
			
			R(jj) = sqrt(  accu( vR.row(jj)%vR.row(jj) ) );	
		}
		
		R6 = square(R)%square(R)%square(R);
		R12 = R6%R6;
		U = accu(4*(1/R12-1/R6));	
		
		dr = delta*(randu(3)-0.5);
		
		r.row(ii) = r.row(ii) + trans(dr);
		
		for (int jj = 0; jj<N && jj != ii; ++jj)
		{
			vR.row(jj) = (r.row(ii) - r.row(jj))  ;
			vR.row(jj) = vR.row(jj) - L*round(vR.row(jj)/L) ;
			
			R(jj) = sqrt(  accu( vR.row(jj)%vR.row(jj) )  );
		
		}
		R6 = square(R)%square(R)%square(R);
		R12 = R6%R6;
		U_ = accu(4*(1/R12-1/R6));
		
		A = std::min(1.,exp((U-U_)/T));
		
		double rand=randu();

		if (rand > A)
		{
			cont = cont +1;
			r.row(ii) = r.row(ii) - trans(dr);	
		}
		else
		{
			cont = cont -1;
		}
		
	}
	return 0;
	
}
