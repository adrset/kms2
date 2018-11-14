#include <cmath>
#include <vector>
#include <ctime>
#include <fstream>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>

#include <TCanvas.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TGraph.h>
using namespace std;

const double a = 0.38;
const double T_0 = 100;
const unsigned int N = 100;
const double k = 0.00831;
const double m = 39.948;
const double R = 0.38;
const double L = 2.3;
const double f = 10000.0;
const double e = 1.0;
const double dx = 1.0 / (double )N;
const double PI = 3.14159265359;
double kappa = 0;
double omega = 0;
double tau   = 0;
double dt = 0.0001;

double pown(double a, int power){
	double ret = a;
	
	for(int ii = 0; ii < power-1; ii++ ){
		ret *= a;
	}

	return ret;
}

Double_t getMax(Double_t tab[], int size){
	Double_t max = tab[0];
	for(int i = 1;i< size ; i++){
		if(tab[i]>max)	
			max = tab[i];
		
	}
	
	return max;
	
}

Double_t getMin(Double_t tab[], int size){
	Double_t min = tab[0];
	for(int i = 1;i< size ; i++){
		if(tab[i]<min)	
			min = tab[i];
		
	}
	
	return min;
	
}

double computeHamiltonian(double hh[], unsigned int k){
	return - 0.5 * ( hh[ k + 1 ] + hh[ k - 1 ] - 2.0 * hh[ k ] ) / (dx * dx) + kappa * ( k * dx - 0.5) * hh[ k ] * sin (omega * tau);
	
}

int main( int argc, char** argv ){
	
	std::ofstream out("output.dat");

	double* re = new double[N + 1];
	double* im = new double[N + 1];
	double* hr = new double[N + 1];
	double* hi = new double[N + 1];
	double* p = new double[N + 1];

	for( unsigned int k = 0 ; k < N + 1; k++ ){
		re[ k ] = sqrt(2) * sin ( PI * k * dx);
		im[ k ] = 0;
		p[ k ] = 0;
	}

	hi [ 0 ] = 0;
	hi [ N ] = 0;
	hr [ 0 ] = 0;
	hr [ N ] = 0;

// Compute Hamiltonians
	for( unsigned int k = 1 ; k < N ; k++ ){

		hi[ k ] = computeHamiltonian(hi, k);

		hr[ k ] = computeHamiltonian(hr, k);

	}

	double tmpRe = 0.;
	for( unsigned int sim = 0; sim < 100000; sim++ ) {

		for( unsigned int k = 1 ; k < N ; k++ ){
			
			// half dt
			re [ k ] = re [ k ] +  hi[ k ] * dt / 2.0;
			hr [ k ] = computeHamiltonian(hr, k);
			
			// imaginary psi and hamiltonian
			im [ k ] = im [ k ] - hr [ k ] * dt;
			hi [ k ] = computeHamiltonian(hi, k);
			
			// real psi and hamiltonian
			re [ k ] = re [ k ] + hi [ k ] * dt / 2.0;
			hr [ k ] = computeHamiltonian(hr, k);
			
			p[ k ]  =  re[ k ] * re[ k ]  + im[ k ] * im[ k ];
		
		}
		
		if(sim % 100 == 0){
			for( unsigned int k = 0 ; k < N + 1 ; k++ ){
				out<< k << "\t" << p[ k ] << std::endl; 
			}
			out<< std::endl<< std::endl;
			
		}


	}



	
	
}























