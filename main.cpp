#include <cmath>
#include <vector>
#include <ctime>
#include <fstream>
#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <cstdlib>

#include <TCanvas.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TGraph.h>
using namespace std;

#include <iostream>
#include <cmath>
#include "Properties.hpp"
//#include <string>

using namespace std;

const double pi = 3.14159265359;

unsigned int N = 0;
double dx = 0;
double dt = 0;
int out_frequency = 0;
int av_frequency = 0;
double kappa = 0;
double omega = 0;
double stat = 0;

double pown(double a, int power)
{
	double ret = a;

	for (int ii = 0; ii < power - 1; ii++)
	{
		ret *= a;
	}

	return ret;
}

Double_t getMax(Double_t tab[], int size)
{
	Double_t max = tab[0];
	for (int i = 1; i < size; i++)
	{
		if (tab[i] > max)
			max = tab[i];
	}

	return max;
}

Double_t getMin(Double_t tab[], int size)
{
	Double_t min = tab[0];
	for (int i = 1; i < size; i++)
	{
		if (tab[i] < min)
			min = tab[i];
	}

	return min;
}

double computeHamiltonian(double hh[], unsigned int k, double tau)
{
	return -0.5 * (hh[k + 1] + hh[k - 1] - 2.0 * hh[k]) / (dx * dx) + kappa * (k * dx - 0.5) * hh[k] * sin(omega * tau);
}

double calculateStationary(unsigned int n, unsigned int k)
{
	return sqrt(2.) * sin(n * pi * k * dx);
}

double calculateNorm(double *im, double *re)
{
	double ret = 1.;
	ret *= dx;
	for (unsigned int ii = 0; ii < N + 1; ii++)
		ret += (re[ii] * re[ii] + im[ii] * im[ii]);
	return ret;
}

double calculateAvPos(double *im, double *re)
{ // calculate average position
	double ret = 1.;
	ret *= dx;
	for (unsigned int ii = 0; ii < N + 1; ii++)
		ret += (ii * dx * (re[ii] * re[ii] + im[ii] * im[ii]));
	return ret;
}

double calculateAvEnergy(double *im, double *re, double *HI, double *HR)
{ // calculate average energy
	double ret = 1.;
	ret *= dx;
	for (unsigned int ii = 0; ii < N + 1; ii++)
		ret += (re[ii] * HR[ii] + im[ii] * HI[ii]);
	return ret * dx;
}

int main(int argc, char *argv[])
{
	if(argc<2)
		return 1;

	std::string file (argv[1]);
	Properties prop("simulation.properties");
	N = prop.getProperty("simulation.N");
	dx = 1. / N;
	stat = prop.getProperty("simulation.n");
	dt = prop.getProperty("simulation.dt");
	out_frequency = prop.getProperty("simulation.out_frequency");
	av_frequency = prop.getProperty("simulation.av_frequency");
	kappa = prop.getProperty("simulation.kappa");
	omega = prop.getProperty("simulation.omega")/2.* pi * pi;

	int steps = prop.getProperty("simulation.steps");
	int frames = prop.getProperty("simulation.frames");
	int whatToDraw = prop.getProperty("simulation.mode");

	std::ofstream out("density.dat");
	double *re = new double[N + 1];
	double *im = new double[N + 1];
	double *p = new double[N + 1];
	double *x = new double[N + 1];
	double *hr = new double[N + 1];
	double *hi = new double[N + 1];

	double **pp = new double*[frames];

	for(int ii=0;ii<frames;ii++){
		pp[ii] = new double[N+1];
		for (unsigned int jj = 0; jj < N + 1; jj++) pp[ii][jj] = 0;
	}
	double tau = 0;

	for (unsigned int ii = 0; ii < N + 1; ii++)
	{
		re[ii] = calculateStationary(stat, ii);
		x[ii] = ii;
		im[ii] = 0;
	}

	for (unsigned int ii = 1; ii < N; ii++)
	{
		hr[ii] = computeHamiltonian(re, ii, tau);
		hi[ii] = computeHamiltonian(im, ii, tau);
	}

	hr[0] = 0.0;
	hr[N] = 0.0;
	hi[0] = 0.0;
	hi[N] = 0.0;
	int currentFrame=0;
	double* energy = new double[steps/av_frequency];
	double* time = new double[steps/av_frequency];
	for (unsigned int sim = 0; sim < steps; sim++)
	{
		//leap frog ?

		for (unsigned int ii = 0; ii < N + 1; ii++)
			re[ii] = re[ii] + hi[ii] * dt / 2.;

		for (unsigned int jj = 1; jj < N; jj++)
			hr[jj] = computeHamiltonian(re, jj, tau);

		for (unsigned int ii = 0; ii < N + 1; ii++)
			im[ii] = im[ii] - hr[ii] * dt;

		for (unsigned int jj = 1; jj < N; jj++)
			hi[jj] = computeHamiltonian(im, jj, tau);

		for (unsigned int ii = 0; ii < N + 1; ii++)
			re[ii] = re[ii] + hi[ii] * dt / 2.;

		for (unsigned int ii = 0; ii < N + 1; ii++)
		{
			p[ii] = re[ii] * re[ii] + im[ii] * im[ii];
		}

		if(sim % av_frequency == 0){
			energy[sim/av_frequency] = calculateAvEnergy(re, im, hr, hi);
			time[sim/av_frequency] = tau;
		}

		if(sim % (steps / frames) == 0){
			std::cout<<"Current frame: "<<currentFrame<<std::endl;
			memcpy(pp[currentFrame++],p , (N+1)*sizeof(double));
		}

		if(sim % out_frequency == 0){
			for (unsigned int ii = 0; ii < N + 1; ii++)
			{
				out<<p[ii]<<std::endl;
			}
			out<<std::endl<<std::endl;
		}
		tau +=dt;
	}

	if(whatToDraw == 1){
		TMultiGraph *mg = new TMultiGraph();
		TString tt = "Density plot; x; ro";

		mg->SetTitle(tt);

		TApplication *rootapp = new TApplication("Pusta", 0, argv);

		TCanvas *tc = new TCanvas("c1", "E", 0, 0, 1400, 400);
		tc->Divide(2,1);
		tc->cd(1);
		TGraph **graphs = new TGraph*[frames]; 
		TLegend* l=new TLegend(0.91,0.425,0.99,0.625);
		l->SetHeader("Density");
		//
		for(int ii=0;ii<frames;ii++){
			graphs[ii] = new TGraph(N + 1, x, pp[ii]);
			TString nam = "graph";
			nam += ii;
			int col = 2 + ii;
			TString nam2 = "";
			nam2 += (ii+1);
			nam2 += "/";
			nam2 += frames;
			nam2+= " * Tau";
			TString color = "l";
			color+= col;
			l->AddEntry(nam,nam2, color);
			graphs[ii]->SetMarkerStyle(19);
			graphs[ii]->SetMarkerColor(col);
			graphs[ii]->SetLineWidth(2);
			graphs[ii]->SetLineColor(col);
			mg->Add(graphs[ii]);
		}
		
		//mg->Add(graph1);
		TGraph *graph2 = new TGraph(steps / av_frequency, time, energy);
		graph2->SetMarkerStyle(19);
		graph2->SetMarkerColor(3);
		graph2->SetLineWidth(2);
		graph2->SetLineColor(3);
		graph2->SetName("graph2");
		graph2->SetTitle("Energy plot; t; E");

		
		mg->Draw("AL");
		l->Draw();
		tc->cd(2);
		graph2->Draw("AL");
		rootapp->Run();
	}else if(whatToDraw == 2){
		system("gnuplot animate.dem");

	}
	for(int ii=0;ii<frames;ii++){
		delete[] pp[ii];
	}

	delete[] pp;
	delete[] hr;
	delete[] hi;
	delete[] re;
	delete[] im;
	delete[] p;
	delete[] x;
}
