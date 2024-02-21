#include "Shallow1D.h"
#include "Shallow2D.h"
#include <algorithm>
#include <stdexcept>
#include <string>
inline double bumpFunc(double xi){
	return 1+0.45*exp(- 5.0*xi*xi);
}
inline double tinyBumpFunc(double xi){
	return 1+0.045*exp(- 20.0*xi*xi);
}
inline double dam2D(double xi, double yi){
	double h = 1.0;
	if (xi>-0.5 && xi<0.5 && yi>-0.5 && yi<0.5 ) h = 2.0;
	return h;
}

int main(int argc, char* argv[]){
	
	/*
	Shallow1D Q1("test1D");
	Q1.setDamBreak(0,3,1);
	Q1.solveSWE();
	*/
	
	
	Shallow2D Q2("test2D");
	Q2.setDomain(-1.0,1.0,-1.0,1.0);
	Q2.setBC(1,1,1,1);
	Q2.setDepthProfile(dam2D);
	Q2.setSimulationTime(0.201);
	Q2.solveSWE();

	return(0);
}

