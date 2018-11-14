#include <cmath>
#include <vector>
#include <ctime>
#include <fstream>
#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>

#include <TCanvas.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TGraph.h>
using namespace std;

const unsigned int N = 100;
const double dx = 1.0 / (double )N;
const double PI = 3.14159265359;
const double hbar = 1.0;
double kappa = 0;
double omega = 0;
double tau   = 0;
double dt = 0.001;
const unsigned int n = 1;

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
		re[ k ] = sqrt(2) * sin ( n * PI * k * dx);
		im[ k ] = 0;
		p[ k ] = 0;
	}

	hi [ 0 ] = 0;
	hi [ N ] = 0;
	hr [ 0 ] = 0;
	hr [ N ] = 0;

// Compute Hamiltonians
	for( unsigned int k = 1 ; k < N ; k++ ){

		hi[ k ] = computeHamiltonian(im, k);

		hr[ k ] = computeHamiltonian(re, k);
	}

	double tmpRe = 0.;
	for( unsigned int sim = 0; sim < 10000; sim++ ) {

		for( unsigned int k = 0 ; k < N+1 ; k++ ){
			
			// half dt
			re [ k ] = re [ k ] +  hi[ k ] * dt / 2.0;
			
		}
		for( unsigned int k = 1 ; k < N ; k++ ){
			hr [ k ] = computeHamiltonian(re, k);
		}
			
			// imaginary psi and hamiltonian
		for( unsigned int k = 0 ; k < N+1 ; k++ ){
			im [ k ] = im [ k ] - hr [ k ] * dt;
		}
		for( unsigned int k = 1 ; k < N ; k++ ){
			hi [ k ] = computeHamiltonian(im, k);
		}
		
			
			// real psi and hamiltonian
		for( unsigned int k = 0 ; k < N + 1; k++ ){
			re [ k ] = re [ k ] + hi [ k ] * dt / 2.0;
		}
			
		for( unsigned int k = 1 ; k < N ; k++ ){
			p[ k ]  =  re[ k ] * re[ k ]  + im[ k ] * im[ k ];
		}
		
		
		
		if(sim % 100 == 0){
			
			//out<< std::endl<< std::endl;
			
		}


	}
	
	for( unsigned int k = 0 ; k < N + 1 ; k++ ){
				out<< re[ k ] << std::endl; 
		}



	
	
}























