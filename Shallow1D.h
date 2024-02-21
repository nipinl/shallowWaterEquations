#ifndef SHALLOW1D_H_
#define SHALLOW1D_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;


class Shallow1D
{	
	double Xmin, Xmax;
	double depth;						//Domain minimum and maximum 
	double gravity;							// acceleration due to gravity
	int Nx;									//Length of vector=No. of control volumes +3
	double cfl;
	int lbc,rbc;							//Boundary condition 0: Outflow or zero gradient, 1: Solid or Reflective bc
	double dt,dx,simTime,endTime;
	vector<double> x;						//x
	vector<double> u1,u2;					//u1=h, u2=hu
	vector<double> u1n,u2n;					// to store old value when updating
	vector<double> u1init,u2init;			// to store initial value for plotting

	ofstream outFile;
	string caseName;

public:
	Shallow1D();
	Shallow1D(string caseName_);
	Shallow1D(double Xmin_, double Xmax_, double depth_, double cells, string caseName_);
	void reMesh(int);
	void details();
	void setCFL(const double cfl_);
	void setSimulationTime(const double endTime_);
	void setGravity(const double gravity_);
	void setDomain(const double Xmin_, const double Xmax_);
	void setDepth(const double depth_);
	void setDepthProfile(double (*depthFunc)(double xi));
	void setDamBreak(const double damLocation, const double upstreamDepth, const double downstreamDepth);
	void setBoundary(const int leftBC,const int rightBC);
	void initialiseDomain(void);
	void initialiseVectors(void);
	void reInitialise(void);
	void solveSWE();
	void applyBC();
	void writeData(const double t,std::ios_base::openmode mode);
	void updateSolution_LaxFrederichs();
	void updateSolution_LF();
	void flux_LaxFrederichs(const int i, double& FL1, double& FL2, double& FR1, double& FR2);
	double getMaxChar(const int i);
	double getGlobalMaxChar(void);
	void printVector(const vector<double>& v);

};

#endif
