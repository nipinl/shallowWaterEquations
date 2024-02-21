#include "Shallow2D.h"

//		 (ghost)        (1)        (2)               (Nx-4)      (Nx-3)      (ghost)
//    |-----------|===========|===========|  ~~~  |===========|===========|-----------|
//    0           1           2           3     Nx-4        Nx-3        Nx-2         Nx-1

Shallow2D::Shallow2D():Nx(103),Ny(103),
									Xmin(-5.0), Xmax(5.0), Ymin(-5.0), Ymax(5.0),
									depth(0.1),
									caseName("defaultCase"),//default no. of control volumes 100. No. of points = mesh+1+2 ghost nodes
									endTime(2.0),
									cfl(0.8),//default value of cfl=1.0
									simTime(0.0),
									lbc(0),
									rbc(0),
                                    bbc(0),
									tbc(0),
									gravity(1.0)

{}
Shallow2D::Shallow2D(string caseName_):Shallow2D()//calling default constructor first and then assign case name
{caseName=caseName_;}
Shallow2D::Shallow2D(double Xmin_, double Xmax_, double Ymin_, double Ymax_,double depth_, double x_cells, double y_cells, string caseName_):Shallow2D(){//calling default constructor first and then assign case name
	Xmin=Xmin_;
	Xmax=Xmax_;
    Ymin=Ymin_;
	Ymax=Ymax_;
	depth=depth_;
	Nx = x_cells+3;
    Ny = y_cells+3;
	caseName=caseName_;
}
//setters
void Shallow2D::reMesh(int x_cells,int y_cells){
	if (x_cells <= 1||y_cells <= 1){
		cout << "At least 1 cell is needed. Setting 1 cell !!" << endl;
        if (x_cells<=1) Nx=x_cells+3;
        if (y_cells<=1) Ny=y_cells+3;
	}
	reInitialise();
}

void Shallow2D::setCFL(double cfl_){cfl = cfl_;}
void Shallow2D::setSimulationTime(double endTime_){
	endTime = endTime_;
	simTime=0;
	}

void Shallow2D::setGravity(double gravity_){gravity = gravity_;}

void Shallow2D::initialiseDomain(){
	simTime=0;
	dx = (Xmax-Xmin)/double(Nx-3);
    dy = (Ymax-Ymin)/double(Ny-3);
	x.push_back(Xmin);					//x[0]: left ghost node		
	for (size_t i = 0; i < Nx-2; i++){	// i: 0 to Nx-3
		x.push_back(Xmin+i*dx); 		//x[1] to x[Nx-2]
	}
	x.push_back(Xmax);					//x[Nx-1]: right ghost node

    y.push_back(Ymin);					//y[0]: bottom ghost node
    for (size_t i = 0; i < Ny-2; i++){	// i: 0 to Ny-3
		y.push_back(Ymin+i*dy); 		//y[1] to y[Ny-2]
	}
    y.push_back(Ymax);					//y[Ny-1]: right ghost node
}
void Shallow2D::initialiseVectors(){
	simTime=0;
    vector<double> h_1d(x.size(), depth);
    vector<double> zero_1d(x.size(), 0.0);	
	for (size_t i = 0; i < y.size(); i++){ 
		u1.push_back(h_1d);
		u2.push_back(zero_1d);
        u3.push_back(zero_1d);
	}
	//inflating the vectors with depth
	u1n		=	u1;
	u1init	=	u1;
	//inflating the vectors with zeros
	u2n		=	u2;
	u2init	=	u2;
    u3n		=	u3; //we can use u3n=u2 and omit inflating u3 in the above loop. but, who cares.
	u3init	=	u3;
}

void Shallow2D::reInitialise(void){
	//clear all vectors
	x.clear();
    y.clear();
	u1.clear();u1n.clear();u1init.clear();	
	u2.clear();u2n.clear();u2init.clear();
    u3.clear();u3n.clear();u3init.clear();
	initialiseDomain();
	initialiseVectors();	
}
void Shallow2D::setDomain(const double Xmin_, const double Xmax_,const double Ymin_, const double Ymax_){
	Xmin=Xmin_;
	Xmax=Xmax_;
    Ymin=Ymin_;
	Ymax=Ymax_;
    if(x.size()== 0){
        initialiseDomain();
        initialiseVectors();
    }
    else{
        reInitialise();
    }
}
void Shallow2D::setDepth(const double depth_){
	depth=depth_;
	reInitialise();
}

void Shallow2D::setDepthProfile(double (*depthFunc)(double xi,double yi)){
	if(x.size()== 0) initialiseDomain();	//check if x,y is initialised
	if(u1.size()== 0) initialiseVectors();	//check if U vectors are initialised
	for (size_t i = 0; i < x.size(); i++){
        for (size_t j = 0; j < y.size(); j++){
		    u1[i][j]	=	depthFunc(x[i],y[j]);
        }
	}
}
/* void Shallow2D::setDamBreak(const double damLocation, const double upstreamDepth, const double downstreamDepth){
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
} */
void Shallow2D::setBC(const int leftBC,const int rightBC, const int bottomBC,const int topBC){
    lbc=leftBC;
    rbc=rightBC;
    bbc=bottomBC;
    tbc=topBC;
}

//getters

double Shallow2D::getMaxCharX(const int i,const int j){
	double h = u1n[i][j];
    double hu = u2n[i][j];
    double u = hu/h;
	double sqgh = sqrt(gravity*h);
	double maxChar = abs(u) + sqgh;
	hu = u2n[i+1][j];
	h = u1n[i+1][j];
	u = hu/h;
	//cout<<"h = " <<h <<" and i+1 is "<<i+1<<endl;
	sqgh = sqrt(gravity*h);
	maxChar = max(maxChar, abs(u) + sqgh);
	//cout<<"maxChar = " <<maxChar<<endl;
	return maxChar;
}
double Shallow2D::getMaxCharY(const int i,const int j){
	double h = u1n[i][j];
    double hv = u3n[i][j];
    double v = hv/h;
	double sqgh = sqrt(gravity*h);
	double maxChar = abs(v) + sqgh;
	hv = u3n[i][j+1];
	h = u1n[i][j+1];
	v = hv/h;
	//cout<<"h = " <<h <<" and j+1 is "<<j+1<<endl;
	sqgh = sqrt(gravity*h);
	maxChar = max(maxChar, abs(v) + sqgh);
	//cout<<"maxChar = " <<maxChar<<endl;
	return maxChar;
}

double Shallow2D::getGlobalMaxCharX(){
	double lamda{0}, speed{0};
	for (size_t i = 1; i < x.size()-1; i++){
        for (size_t j = 1; j < y.size()-1; j++){
            speed = getMaxCharX(i,j);
		    if (lamda < speed) lamda = speed;
        }
	}
	return lamda;
}
double Shallow2D::getGlobalMaxCharY(){
	double lamda{0}, speed{0};
	for (size_t i = 1; i < x.size()-1; i++){
        for (size_t j = 1; j < y.size()-1; j++){
            speed = getMaxCharY(i,j);
		    if (lamda < speed) lamda = speed;
        }
	}
	return lamda;
}

void Shallow2D::solveSWE(){

	if(x.size()== 0) initialiseDomain();
	if(u1.size()== 0) initialiseVectors();
	details();

	bool milesToGo=true;
	int iterations=1;
	writeData(simTime,std::ofstream::trunc);
    writeMidYData(simTime,std::ofstream::trunc);
	
	
	applyBC();
	while(milesToGo==true){
	//while(iterations<2){
		iterations++;

		// Calculate time step
		dt = 0.5*cfl*min(dx/getGlobalMaxCharX(),dy/getGlobalMaxCharY());
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
        u3n = u3;						
		
		updateSolution_LaxFrederichs();
		
		applyBC();

        writeMidYData(simTime, ios::app);
        writeMidXData(simTime, ios::app);
		writeData(simTime, ios::app);
	}
    //for plotting midlplane
	string iterFileName = caseName+".r";
	outFile.open(iterFileName,ios::app);
	outFile <<(Nx-2)*(Ny-2)<<"\t"<<Nx-2<<"\t"<<Ny-2<<"\t"<<iterations<<"\t"<<Xmin<<"\t"<<Xmax<<"\t"<<Ymin<<"\t"<<Ymax<<endl;
	outFile.close();
}
void Shallow2D::applyBC(){
	if (lbc==0){//Zero-gradient
        for (size_t j = 0; j < y.size(); j++)
        {
            u1[0][j] =	u1[1][j];
			u2[0][j] =	u2[1][j];
            u3[0][j] =	u3[1][j];
        }	
	}
	if(lbc==1){//solid wall or reflective bc
        for (size_t j = 0; j < y.size(); j++)
        {
            u1[0][j] =	u1[1][j];
			u2[0][j] =	-u2[1][j];
            u3[0][j] =	u3[1][j];
        }

	}
	if (rbc==0){//Zero-gradient
        for (size_t j = 0; j < y.size(); j++)
        {
            u1[Nx-1][j] =	u1[Nx-2][j];
			u2[Nx-1][j] =	u2[Nx-2][j];
            u3[Nx-1][j] =	u3[Nx-2][j];
        }
	}
	if(rbc==1){//solid wall or reflective bc
        for (size_t j = 0; j < y.size(); j++)
        {
            u1[Nx-1][j] =	u1[Nx-2][j];
			u2[Nx-1][j] =	-u2[Nx-2][j];
            u3[Nx-1][j] =	u3[Nx-2][j];
        }
	}
//-------------------------------------------------------

    if (bbc==0){//Zero-gradient
        for (size_t i = 0; i < x.size(); i++)
        {
            u1[i][0] =	u1[i][1];
			u2[i][0] =	u2[i][1];
            u3[i][0] =	u3[i][1];
        }	
	}
	if(bbc==1){//solid wall or reflective bc
        for (size_t i = 0; i < x.size(); i++)
        {
            u1[i][0] =	u1[i][1];
			u2[i][0] =	u2[i][1];
            u3[i][0] =	-u3[i][1];
        }

	}
	if (tbc==0){//Zero-gradient
        for (size_t i = 0; i < x.size(); i++)
        {
            u1[i][Ny-1] =	u1[i][Ny-2];
			u2[i][Ny-1] =	u2[i][Ny-2];
            u3[i][Ny-1] =	u3[i][Ny-2];
        }
	}
	if(tbc==1){//solid wall or reflective bc
        for (size_t i = 0; i < x.size(); i++)
        {
            u1[i][Ny-1] =	u1[i][Ny-2];
			u2[i][Ny-1] =	u2[i][Ny-2];
            u3[i][Ny-1] =	-u3[i][Ny-2];
        }
	}	
}
void Shallow2D::updateSolution_LaxFrederichs(){
	double FL1{0},FL2{0},FL3{0},FR1{0},FR2{0},FR3{0};
    double GL1{0},GL2{0},GL3{0},GR1{0},GR2{0},GR3{0};
	for(size_t i=1;i<Nx-1;i++){
        for(size_t j=1;j<Ny-1;j++){	
		//Get left and right fluxes (Lax-Friedrichs)
            flux_LaxFrederichs(i,j, FL1,FL2,FL3,FR1,FR2,FR3,GL1,GL2,GL3,GR1,GR2,GR3);
            u1[i][j] = u1n[i][j] + (dt/dx)* ( FL1 - FR1 ) + (dt/dy)* ( GL1 - GR1 );
            u2[i][j] = u2n[i][j] + (dt/dx)* ( FL2 - FR2 ) + (dt/dy)* ( GL2 - GR2 );
            u3[i][j] = u3n[i][j] + (dt/dx)* ( FL3 - FR3 ) + (dt/dy)* ( GL3 - GR3 );
	    }
    }
}

void Shallow2D::flux_LaxFrederichs(const int i, const int j, double& FL1, double& FL2, double& FL3, double& FR1, double& FR2, double& FR3
							, double& GL1, double& GL2, double& GL3, double& GR1, double& GR2, double& GR3){
	//F=[f1, f2]', U = [u1, u2]'
	//l:i, r:i+1
	double q1 	= 	u1n[i][j];
	double q1p1	= 	u1n[i+1][j];
	double q1m1	= 	u1n[i-1][j];
	double q2 	= 	u2n[i][j];
	double q2p1	= 	u2n[i+1][j];
	double q2m1	= 	u2n[i-1][j];
    double q3 	= 	u3n[i][j];
	double q3p1	= 	u3n[i+1][j];
	double q3m1	= 	u3n[i-1][j];
	
	double f1 	= 	q2;
	double f1p1	= 	q2p1;
	double f1m1	= 	q2m1;
	double f2 	= 	q2*q2/q1+0.5*gravity*q1*q1;
	double f2p1	= 	q2p1*q2p1/q1p1+0.5*gravity*q1p1*q1p1;
	double f2m1	= 	q2m1*q2m1/q1m1+0.5*gravity*q1m1*q1m1;
    double f3   =   q2*q3/q1;
    double f3p1 =   q2p1*q3p1/q1p1;
    double f3m1 =   q2m1*q3m1/q1m1;

	double maxCharLeft = getMaxCharX(i-1,j);
	double maxCharRight = getMaxCharX(i,j);

	FL1 = 0.5*(f1m1+f1) - 0.5*maxCharLeft * (q1 - q1m1);
	FL2 = 0.5*(f2m1+f2) - 0.5*maxCharLeft * (q2 - q2m1);
    FL3 = 0.5*(f3m1+f3) - 0.5*maxCharLeft * (q3 - q3m1);

	FR1 = 0.5*(f1+f1p1) - 0.5*maxCharRight * (q1p1 - q1);
	FR2 = 0.5*(f2+f2p1) - 0.5*maxCharRight * (q2p1 - q2);
    FR3 = 0.5*(f3+f3p1) - 0.5*maxCharRight * (q3p1 - q3);
//y direction. q1,q2,q3 are same hence only *p1, *m1 are overwritten
	q1p1	= 	u1n[i][j+1];
	q1m1	= 	u1n[i][j-1];
	q2p1	= 	u2n[i][j+1];
	q2m1	= 	u2n[i][j-1];
	q3p1	= 	u3n[i][j+1];
	q3m1	= 	u3n[i][j-1];

    double g1 	= 	q3;
	double g1p1	= 	q3p1;
	double g1m1	= 	q3m1;
	double g2 	= 	q2*q3/q1;
	double g2p1	= 	q2p1*q3p1/q1p1;
	double g2m1	= 	q2m1*q3m1/q1m1;
    double g3   =   q3*q3/q1+0.5*gravity*q1*q1;
    double g3p1 =   q3p1*q3p1/q1p1+0.5*gravity*q1p1*q1p1;
    double g3m1 =   q3m1*q3m1/q1m1+0.5*gravity*q1m1*q1m1;

    maxCharLeft = getMaxCharY(i,j-1);
	maxCharRight = getMaxCharY(i,j);

	GL1 = 0.5*(g1m1+g1) - 0.5*maxCharLeft * (q1 - q1m1);
	GL2 = 0.5*(g2m1+g2) - 0.5*maxCharLeft * (q2 - q2m1);
    GL3 = 0.5*(g3m1+g3) - 0.5*maxCharLeft * (q3 - q3m1);

	GR1 = 0.5*(g1+g1p1) - 0.5*maxCharRight * (q1p1 - q1);
	GR2 = 0.5*(g2+g2p1) - 0.5*maxCharRight * (q2p1 - q2);
    GR3 = 0.5*(g3+g3p1) - 0.5*maxCharRight * (q3p1 - q3);
}


void Shallow2D::writeData(double t,std::ios_base::openmode mode){
	string fileName=caseName+".u2";
	
	outFile.open(fileName,mode);
		for(size_t i=1;i<Nx-1;i++){
            for(size_t j=1;j<Ny-1;j++){
				outFile <<t<<","<< x[i]<<","<< y[j]<<","<<u1[i][j]<<","<< u2[i][j]<<","<< u3[i][j]
                <<","<< u2[i][j]/u1[i][j]<<","<< u3[i][j]/u1[i][j]<<endl;
		    }
        }
	outFile.close();
}
void Shallow2D::writeMidYData(double t,std::ios_base::openmode mode){
	string fileName=caseName+".ux";
	
	outFile.open(fileName,mode);
        int j = ceil(0.5*(Ny-2));
		for(size_t i=1;i<Nx-1;i++){
				outFile <<t<<","<< x[i]<<","<< u1init[i][j]<<","<< u2init[i][j]<<","<< u2init[i][j]/u1init[i][j]
						<<","<< u1[i][j]<<","<< u2[i][j]<<","<< u2[i][j]/u1[i][j]<<endl;
		}
	outFile.close();
}
void Shallow2D::writeMidXData(double t,std::ios_base::openmode mode){
	string fileName=caseName+".uy";
	
	outFile.open(fileName,mode);
        int i = ceil(0.5*(Nx-2));
		for(size_t j=1;j<Ny-1;j++){
				outFile <<t<<","<< y[j]<<","<< u1init[i][j]<<","<< u2init[i][j]<<","<< u2init[i][j]/u1init[i][j]
						<<","<< u1[i][j]<<","<< u2[i][j]<<","<< u2[i][j]/u1[i][j]<<endl;
		}
	outFile.close();
}
void Shallow2D::details(){
	cout<<"Problem: "<<caseName<<endl;
	cout << "Total no. of cells  = "<<(Nx-3)*(Ny-3) <<". "<<Nx-3<<" in x and "<<Ny-3<<" in y"<<endl;
	cout << "CFL = "<<cfl <<" dx = "<<dx<<" dy = "<<dy<< endl;
	cout<<"Numerical flux is Lax Frederick Local";
	cout<<endl;
} 
void Shallow2D::printVector(const vector<double>& v){
	cout<<endl;
	for (size_t i = 0; i < x.size(); i++){
		cout<<v[i]<<" ";
	}
	cout<<endl;
}
void print2dVector(vector<vector<double>> const &v)
{
	for (auto i : v)
	{
		//i is now an 1D vector
		for (auto j : i)
		{
			cout << j << " ";
		}
		cout << endl;
	}
}
