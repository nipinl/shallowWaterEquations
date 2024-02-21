#ifndef Shallow2D_H_
#define Shallow2D_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;


class Shallow2D
{	
	double Xmin, Xmax, Ymin, Ymax;
	double depth;						        //Domain minimum and maximum 
	double gravity;							    // acceleration due to gravity
	int Nx,Ny;									//Length of vector=No. of control volumes +3
	double cfl;
	int lbc,rbc,tbc,bbc;						//Boundary condition 0: Outflow or zero gradient, 1: Solid or Reflective bc
	double dt,dx,dy,simTime,endTime;
	vector<double> x,y;						    //x
	vector<vector<double>> u1,u2,u3;					    //u1=h, u2=hu
	vector<vector<double>> u1n,u2n,u3n;					    // to store old value when updating
	vector<vector<double>> u1init,u2init,u3init;			    // to store initial value for plotting

	ofstream outFile;
	string caseName;

public:
	Shallow2D();
	Shallow2D(string caseName_);
	Shallow2D(double Xmin_, double Xmax_, double Ymin_, double Ymax_, double depth_, double x_cells, double y_cells, string caseName_);
	void reMesh(int x_cells,int y_cells);
	void details();
	void setCFL(const double cfl_);
	void setSimulationTime(const double endTime_);
	void setGravity(const double gravity_);
	void setDomain(const double Xmin_, const double Xmax_,const double Ymin_, const double Ymax_);
	void setDepth(const double depth_);
	void setDepthProfile(double (*depthFunc)(double xi,double yi));
    void setBC(const int leftBC,const int rightBC, const int bottomBC,const int topBC);
	void initialiseDomain(void);
	void initialiseVectors(void);
	void reInitialise(void);
	void solveSWE();
	void applyBC();
	void writeData(const double t,std::ios_base::openmode mode);
	void writeMidXData(double t,std::ios_base::openmode mode);
	void writeMidYData(double t,std::ios_base::openmode mode);
	void updateSolution_LaxFrederichs();
	void updateSolution_LF();
	void flux_LaxFrederichs(const int i, const int j, double& FL1, double& FL2, double& FL3, double& FR1, double& FR2, double& FR3
							, double& GL1, double& GL2, double& GL3, double& GR1, double& GR2, double& GR3);
	double getMaxCharX(const int i,const int j);
	double getMaxCharY(const int i,const int j);
	double getGlobalMaxCharX(void);
	double getGlobalMaxCharY(void);
	void printVector(const vector<double>& v);
};
void print2dVector(vector<vector<double>> const &v);
#endif
