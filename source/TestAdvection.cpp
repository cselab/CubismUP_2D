//
//  TestAdvection.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestAdvection.h"
#include "ProcessOperatorsOMP.h"
#include "OperatorVorticity.h"
#include <sstream>
#include <cmath>

#include "CoordinatorAdvection.h"

/*
 class BS4
 {
 public:
	static inline Real eval(Real x)
	{
 const Real t = fabs(x);
 
 if (t>2) return 0;
 
 if (t>1) return pow(2-t,3)/6;
 
 return (1 + 3*(1-t)*(1 + (1-t)*(1 - (1-t))))/6;
	}
 };
 */

void TestAdvection::_icLinear()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				double p[3];
				info.pos(p, ix, iy);
				
				double rho = 0;
				if (ix==FluidBlock::sizeX/2 && iy==FluidBlock::sizeY/2)
					rho = 1;
				b(ix, iy).rho = rho;//abs(p[0]-.5);
				b(ix, iy).u   = 1;
				b(ix, iy).v   = 1;
				b(ix, iy).chi = 0;
			}
	}
	
	
	stringstream ss;
	ss << path2file << "-IC.vti" ;
	//cout << ss.str() << endl;
	
	dumper.Write(*grid, ss.str());
}

void TestAdvection::_icVortex()
{
	const double center[2] = {.5,.5};
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				double p[3];
				info.pos(p, ix, iy);
				
				p[0] = p[0]*2.-1.;
				p[1] = p[1]*2.-1.;
				
				const Real r = sqrt(p[0]*p[0] + p[1]*p[1]);
				const Real invR = 1/r;
				
				
				b(ix, iy).rho = r;
				b(ix, iy).u   = -p[1];// sin(p[1])*cos(r*M_PI/2)*invR;
				b(ix, iy).v   =  p[0];//-sin(p[0])*cos(r*M_PI/2)*invR;
				b(ix, iy).chi = 0;
				/*
				if (r>.5)
				{
					b(ix,iy).rho = 1;
				}
				/*
				 p[0] = p[0]*2.-1.;
				 p[1] = p[1]*2.-1.;
				 
				 const Real IrI  = sqrt(p[0]*p[0] + p[1]*p[1]);
				 const double strength = 100./(1+IrI*IrI)*BS4::eval(IrI/0.5*2.5);
				 
				 b(ix, iy).rho = p[0];
				 b(ix, iy).u   = -p[1]*strength;
				 b(ix, iy).v   =  p[0]*strength;
				 b(ix, iy).chi = 0;
				 
				/*/
				/*
				const Real dx = p[0] - center[0];
				const Real dy = p[1] - center[1];
				const Real dist = sqrt(dx*dx + dy*dy);
				
				b(ix, iy).rho = abs(p[0]-.5);
				if (dist <= .4)
				{
					const Real amplitude = .5*(cos(dist*5*M_PI+M_PI)+1);
					b(ix, iy).u = -amplitude*dy/dist;
					b(ix, iy).v =  amplitude*dx/dist;
				}
				else
				{
					b(ix, iy).u = 0;
					b(ix, iy).v = 0;
				}
				b(ix, iy).chi = 0;
				
				b(ix, iy).tmpU = 0;
				b(ix, iy).tmpV = 0;
				b(ix, iy).tmp  = 0;
				
				//*/
			}
	}
	
	
	stringstream ss;
	ss << path2file << "-IC.vti";
	
	dumper.Write(*grid, ss.str());
	
	
	const int sizeX = bpd * FluidBlock::sizeX;
	const int sizeY = bpd * FluidBlock::sizeY;
	vorticityIC = new Layer(sizeX,sizeY,1);
	processOMP<Lab, OperatorVorticity>(*vorticityIC,vInfo,*grid);
	stringstream sVort;
	sVort << path2file << "Vorticity-IC.vti";
	dumpLayer2VTK(0,sVort.str(),*vorticityIC,1);
}

void TestAdvection::_icBurger()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				double p[3];
				info.pos(p, ix, iy);
				
				b(ix, iy).rho = 1;
				b(ix, iy).u   = ix<FluidBlock::sizeX/2 ? p[0] : (1-p[0]);//1-p[0];//1-cos(p[0]*M_PI*2);
				b(ix, iy).v   = 0;
				b(ix, iy).chi = 0;
			}
	}
	
	
	stringstream ss;
	ss << path2file << "-IC.vti" ;
	//cout << ss.str() << endl;
	
	dumper.Write(*grid, ss.str());
}

TestAdvection::TestAdvection(const int argc, const char ** argv, int testCase, const int bpd, const double dt, const double tEnd) : Test(argc, argv), time(0), testCase(testCase), bpd(bpd), dt(dt), tEnd(tEnd)
{
	grid = new FluidGrid(bpd,bpd,1);
	
	// setup initial condition
	if (testCase==0)
	{
		// output settings
		path2file = parser("-file").asString("../data/testAdvectionLinear");
		_icLinear();
	}
	else if (testCase==1)
	{
		// output settings
		path2file = parser("-file").asString("../data/testAdvectionVortex");
		_icVortex();
		//path2file = parser("-file").asString("../data/testAdvectionBurger");
		//_icBurger();
	}
	else
	{
		cout << "unknown test case - aborting!\n";
		abort();
	}
}

TestAdvection::~TestAdvection()
{
	delete grid;
}

void TestAdvection::run()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	/*
	 if (testCase==0)
		cout << "Using dt " << dt << " (CFL time step: " << vInfo[0].h_gridpoint/1. << ")\n";
	 else
		cout << "Using dt " << dt << " (CFL time step: " << vInfo[0].h_gridpoint/1. << ")\n";
	 //*/
	
	int step;
	const int nsteps = tEnd/dt;
	
	//CoordinatorAdvection<Lab> coordAdvection(grid);
	CoordinatorTransport<Lab> coordTransport(grid);
	
	while (time<=tEnd)
	{
		//coordAdvection(dt);
		coordTransport(dt);
		
		//dump some time steps every now and then
		/*
		if(step % 10 == 0)
		{
			stringstream ss;
			ss << path2file << "-" << step << ".vti" ;
			
			dumper.Write(*grid, ss.str());
		}
		 */
		
		time += dt;
		step++;
	}
	
	stringstream ss;
	ss << path2file << "-test" << testCase << "-bpd" << bpd << ".vti";
	dumper.Write(*grid, ss.str());
}

void TestAdvection::check()
{
	const double center[2] = {.5,.5};
	
	cout << "\tErrors (uLinf, uL1, uL2):\t";
	double uLinf = 0.;
	double uL1 = 0.;
	double uL2 = 0.;
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const double dh = vInfo[0].h_gridpoint;
	
	
	const int sizeX = bpd * FluidBlock::sizeX;
	const int sizeY = bpd * FluidBlock::sizeY;
	
	if (testCase==0)
	{
#pragma omp parallel for reduction(max:uLinf) reduction(+:uL1) reduction(+:uL2)
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					double p[3];
					info.pos(p, ix, iy);
					
					double uError, vError;
					uError = b(ix, iy).u - 1;
					vError = b(ix, iy).v - 1;
					
					uLinf = max(uLinf,abs(uError));
					uL1 += abs(uError);
					uL2 += uError*uError;
				}
		}
	}
	else
	{
#pragma omp parallel for reduction(max:uLinf) reduction(+:uL1) reduction(+:uL2)
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					double p[3];
					info.pos(p, ix, iy);
					
					p[0] = p[0]*2.-1.;
					p[1] = p[1]*2.-1.;
					
					Real r = sqrt(p[0]*p[0] + p[1]*p[1]);
					
					double error = b(ix, iy).rho - r;
					b(ix,iy).chi = error;
					
					uLinf = max(uLinf,abs(error));
					uL1 += abs(error);
					uL2 += error*error;
				}
		}
		/*
		 #pragma omp parallel for reduction(max:uLinf) reduction(+:uL1) reduction(+:uL2)
		 for (int iy=0; iy<sizeY; iy++)
			for (int ix=0; ix<sizeX; ix++)
			{
		 double error = vorticity(ix,iy)-(*vorticityIC)(ix,iy);
		 vorticityDiff(ix,iy) = error;
		 
		 uLinf = max(uLinf,abs(error));
		 uL1 += abs(error);
		 uL2 += error*error;
			}
		 */
		/*
#pragma omp parallel for reduction(max:uLinf) reduction(+:uL1) reduction(+:uL2)
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					// 1D Burger's
					double p[3];
					info.pos(p, ix, iy);
					double pIC = p[0] - b(ix,iy).u * time; // why this factor 2?
					double error = b(ix, iy).u - (1-cos(pIC*M_PI*2));
					
					b(ix,iy).u = 1-cos(pIC*M_PI*2);
					
					uLinf = max(uLinf,abs(error));
					uL1 += abs(error);
					uL2 += error*error;
				}
		}
		 */
	}
	
	//stringstream sVort;
	//sVort << path2file << "VorticityDiff-" << bpd << ".vti";
	//dumpLayer2VTK(0,sVort.str(),vorticityDiff,1);
	
	stringstream ssol;
	ssol << path2file << "-solution" << testCase << "-bpd" << bpd << ".vti";
	dumper.Write(*grid, ssol.str());
	
	uL1 *= dh*dh;
	uL2 = sqrt(uL2)*dh;
	cout << "\t" << uLinf << "\t" << uL1 << "\t" << uL2 << endl;
	myfile << sizeX << " " << uLinf << " " << uL1 << " " << uL2 << endl;
}
