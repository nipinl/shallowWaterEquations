#include "Shallow1D.h"

//		 (ghost)        (1)        (2)               (Nx-4)      (Nx-3)      (ghost)
//    |-----------|===========|===========|  ~~~  |===========|===========|-----------|
//    0           1           2           3     Nx-4        Nx-3        Nx-2         Nx-1

Shallow1D::Shallow1D():Nx(103),
									Xmin(-5.0), Xmax(5.0),
									depth(0.1),
									caseName("defaultCase"),//default no. of control volumes 100. No. of points = mesh+1+2 ghost nodes
									endTime(2.0),
									cfl(0.8),//default value of cfl=1.0
									simTime(0.0),
									lbc(0),
									rbc(0),
									gravity(1.0)

{}
Shallow1D::Shallow1D(string caseName_):Shallow1D()//calling default constructor first and then assign case name
{caseName=caseName_;}
Shallow1D::Shallow1D(double Xmin_, double Xmax_, double depth_, double cells, string caseName_):Shallow1D(){//calling default constructor first and then assign case name
	Xmin=Xmin_;
	Xmax=Xmax_;
	depth=depth_;
	Nx = cells+3;
	caseName=caseName_;
}
//setters
void Shallow1D::reMesh(int cells){
	if (cells <= 1){
		cout << "At least 1 cell is needed. Setting 1 cell !!" << endl;
	}
	Nx=cells+3;
	reInitialise();
}

void Shallow1D::setCFL(double cfl_){cfl = cfl_;}
void Shallow1D::setSimulationTime(double endTime_){
	endTime = endTime_;
	simTime=0;
	}

void Shallow1D::setGravity(double gravity_){gravity = gravity_;}

void Shallow1D::initialiseDomain(){
	simTime=0;
	dx = (Xmax-Xmin)/double(Nx-3);
	x.push_back(Xmin);					//x[0]: left ghost node	
	for (size_t i = 0; i < Nx-2; i++){	// i: 0 to Nx-3
		x.push_back(Xmin+i*dx); 		//x[1] to x[Nx-2]
	}
	x.push_back(Xmax);					//x[Nx-1]: right ghost node
}
void Shallow1D::initialiseVectors(){
	simTime=0;	
	for (size_t i = 0; i < x.size(); i++){ 
		u1.push_back(depth);
		u2.push_back(0);
	}
	//inflating the vectors with depth
	u1n		=	u1;
	u1init	=	u1;
	//inflating the vectors with zeros
	u2n		=	u2;
	u2init	=	u2;
}

void Shallow1D::reInitialise(void){
	//clear all vectors
	x.clear();
	u1.clear();u1n.clear();u1init.clear();	
	u2.clear();u2n.clear();u2init.clear();
	initialiseDomain();
	initialiseVectors();	
}
void Shallow1D::setDomain(const double Xmin_, const double Xmax_){
	Xmin=Xmin_;
	Xmax=Xmax_;
	initialiseDomain();
}
void Shallow1D::setDepth(const double depth_){
	depth=depth_;
	reInitialise();
}
void Shallow1D::setDepthProfile(double (*depthFunc)(double xi)){
	if(x.size()== 0) initialiseDomain();	//check if x is initialised
	if(u1.size()== 0) initialiseVectors();	//check if u vectors are initialised
	for (size_t i = 0; i < x.size(); i++){
		u1[i]	=	depthFunc(x[i]);
	}
}
void Shallow1D::setDamBreak(const double damLocation, const double upstreamDepth, const double downstreamDepth){
	double xd = damLocation;
	if(x.size()== 0) initialiseDomain();	//check if x is initialised
	if(u1.size()== 0) initialiseVectors();	//check if u vectors are initialised
	if (damLocation<=Xmin || damLocation>=Xmax) xd = Xmin+0.3*(Xmax-Xmin);
	for (size_t i = 0; i < x.size(); i++){
		if (x[i]<=xd){
			u1[i]	=	upstreamDepth;
		}
		else{
			u1[i]	=	downstreamDepth;
		}		
	}
	//copy u1 to u1n and u1init
	u1n		=	u1;
	u1init	=	u1;
	//velocity is zero for damBreak. Hence purging u2, u2n, u2init to zero.
	std::fill(u2.begin(), u2.end(), 0.0);
	std::fill(u2n.begin(), u2n.end(), 0.0);
	std::fill(u2init.begin(), u2init.end(), 0.0);
}
void Shallow1D::setBoundary(const int leftBC,const int rightBC){
    lbc=leftBC;
    rbc=rightBC;
}

//getters
	
double Shallow1D::getMaxChar(const int i){
	double hu = u2n[i];
	double h = u1n[i];
	double u = hu/h;
	double sqgh = sqrt(gravity*h);
	double maxChar = max(abs(u + sqgh) , abs(u - sqgh));
	hu = u2n[i+1];
	h = u1n[i+1];
	u = hu/h;
	//cout<<"h = " <<h <<" and i+1 is "<<i+1<<endl;
	sqgh = sqrt(gravity*h);
	maxChar = max(maxChar, abs(u + sqgh));
	maxChar = max(maxChar, abs(u - sqgh));
	//cout<<"maxChar = " <<maxChar<<endl;

	return maxChar;
}

double Shallow1D::getGlobalMaxChar(){
	double lamda{0};
	for (size_t i = 1; i < x.size()-1; i++){
		if (lamda < getMaxChar(i)) lamda = getMaxChar(i);
	}
	return lamda;
}

void Shallow1D::solveSWE(){

	if(x.size()== 0) initialiseDomain();
	if(u1.size()== 0) initialiseVectors();
	details();

	bool milesToGo=true;
	int iterations=1;
	writeData(simTime,std::ofstream::trunc);
	
	
	applyBC();
	while(milesToGo==true){
	//while(iterations<2){
		iterations++;

		// Calculate time step
		double maxChar = getGlobalMaxChar();
		dt = cfl*dx/maxChar;
		if (simTime+dt>=endTime){
			dt=min(endTime-simTime,dt);
			milesToGo=false;
		}
		
		//update simulation time
		simTime+=dt;
		cout<<simTime<<" s"<<endl;

		//copy variables before changing
		u1n = u1;
		u2n = u2;						
		
		//updateSolution_LF();
		/* cout<<"u1"<<endl;
		printVector(u1);
		cout<<"u2"<<endl;
		printVector(u2); */
		updateSolution_LaxFrederichs();
		
		applyBC();


		writeData(simTime, ios::app);
	}
	string iterFileName = caseName+".r";
	outFile.open(iterFileName,ios::app);
	outFile <<Nx-2<<"\t"<<iterations<<"\t"<<Xmin<<"\t"<<Xmax<<"\t"<<endl;
	outFile.close();
}
void Shallow1D::applyBC(){
	if (lbc==0){//Zero-gradient
			u1[0] =	u1[1];
			u2[0] =	u2[1];
		}
		else if(lbc==1){//solid wall or reflective bc
			u1[0] =	u1[1];
			u2[0] =	-u2[1];
		}
		if (rbc==0){//Zero-gradient
			u1[Nx-1] =	u1[Nx-2];
			u2[Nx-1] =	u2[Nx-2];
		}
		else if(rbc==1){//solid wall or reflective bc
			u1[Nx-1] =	u1[Nx-2];
			u2[Nx-1] =	-u2[Nx-2];
		}	
}
void Shallow1D::updateSolution_LaxFrederichs(){
	double FL1{0},FL2{0},FR1{0},FR2{0};
	for(size_t i=1;i<Nx-1;i++){	
		//Get left and right fluxes (Lax-Friedrichs)
		flux_LaxFrederichs(i, FL1,FL2,FR1,FR2);
		u1[i] = u1n[i] + (dt/dx)* ( FL1 - FR1 );
		u2[i] = u2n[i] + (dt/dx)* ( FL2 - FR2 );
	}
}
void Shallow1D::updateSolution_LF(){
/* 	double F1{0},F2{0};

//	un1 = h;                 % Function to calculate initial
//	un2 = h.*v;              % dependent variables.

//	F1 = h.*v;
//	F2 = h.*v.^2 + (g.*h.^2)/2;

	for(size_t i=1;i<Nx-1;i++){	
		//Get left and right fluxes (Lax-Friedrichs)
		F1 = (u2n[i+1]*u2n[i+1]/u1n[i+1] + 0.5*gravity*u1n[i+1]*u1n[i+1]);
		F2 = (u2n[i-1]*u2n[i-1]/u1n[i-1] + 0.5*gravity*u1n[i-1]*u1n[i-1]);
		u1[i] = 0.5*(u1n[i-1]+u1n[i+1]) - 0.5* (dt/dx)* ( u2n[i+1] - u2n[i-1] );
		u2[i] = 0.5*(u2n[i-1]+u2n[i+1]) - 0.5* (dt/dx)* ( F1 - F2 );
	} */
}
void Shallow1D::flux_LaxFrederichs(const int i, double& FL1, double& FL2, double& FR1, double& FR2){
	//F=[f1, f2]', U = [u1, u2]'
	//l:i, r:i+1
	double q1 	= 	u1n[i];
	double q1p1	= 	u1n[i+1];
	double q1m1	= 	u1n[i-1];
	double q2 	= 	u2n[i];
	double q2p1	= 	u2n[i+1];
	double q2m1	= 	u2n[i-1];
	
	double f1 	= 	q2;
	double f1p1	= 	q2p1;
	double f1m1	= 	q2m1;
	double f2 	= 	q2*q2/q1+0.5*gravity*q1*q1;
	double f2p1	= 	q2p1*q2p1/q1p1+0.5*gravity*q1p1*q1p1;
	double f2m1	= 	q2m1*q2m1/q1m1+0.5*gravity*q1m1*q1m1;

	double maxCharLeft = getMaxChar(i-1);
	double maxCharRight = getMaxChar(i);

	FL1 = 0.5*(f1m1+f1) - 0.5*maxCharLeft * (q1 - q1m1);
	FL2 = 0.5*(f2m1+f2) - 0.5*maxCharLeft * (q2 - q2m1);

	FR1 = 0.5*(f1+f1p1) - 0.5*maxCharRight * (q1p1 - q1);
	FR2 = 0.5*(f2+f2p1) - 0.5*maxCharRight * (q2p1 - q2);
}


void Shallow1D::writeData(double t,std::ios_base::openmode mode){
	string fileName=caseName+".u";
	
	outFile.open(fileName,mode);
		for(size_t i=1;i<Nx-1;i++){
				outFile <<t<<","<< x[i]<<","<< u1init[i]<<","<< u2init[i]<<","<< u2init[i]/u1init[i]
						<<","<< u1[i]<<","<< u2[i]<<","<< u2[i]/u1[i]<<endl;
		}
	outFile.close();
}
void Shallow1D::details(){
	cout<<"Problem: "<<caseName<<endl;
	cout << "No. of cells  = "<<Nx-3 << endl;
	cout << "CFL = "<<cfl <<" dx = "<<dx<<" dt = "<<dt<< endl;
	cout<<"Numerical flux is Lax Frederick Local";
	cout<<endl;
} 
void Shallow1D::printVector(const vector<double>& v){
	cout<<endl;
	for (size_t i = 0; i < x.size(); i++){
		cout<<v[i]<<" ";
	}
	cout<<endl;
}
