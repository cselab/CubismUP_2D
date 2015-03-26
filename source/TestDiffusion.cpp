//
//  TestDiffusion.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestDiffusion.h"
#include "ProcessOperatorsOMP.h"
#include <sstream>
#include <cmath>

double TestDiffusion::_analytical(double px, double py, double t)
{
	return sin(px*2.*M_PI) * sin(py*2.*M_PI) * exp(-8.*nu*M_PI*M_PI*t);
}

void TestDiffusion::_ic()
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
				
				b(ix, iy).rho = 0;
				b(ix, iy).u   = _analytical(p[0],p[1],0);
				b(ix, iy).v   = 0;
				b(ix, iy).chi = 0;
			}
	}
	
	
	stringstream ss;
	ss << path2file << "-IC.vti" ;
	//cout << ss.str() << endl;
	
	dumper.Write(*grid, ss.str());
}

TestDiffusion::TestDiffusion(const int argc, const char ** argv, const int bpd) : Test(argc, argv), bpd(bpd), time(0)
{
	// test settings
	nu = parser("-nu").asDouble(1);
	
	// output settings
	path2file = parser("-file").asString("../data/testDiffusion");
	
	grid = new FluidGrid(bpd,bpd,1);
	
	// setup initial condition
	_ic();
}

TestDiffusion::~TestDiffusion()
{
	delete grid;
}

void TestDiffusion::run()
{
	time = 0;
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	// TODO: choose dt (Fourier)
	const double dt = 1e-7;
	//cout << "Using dt " << dt << " (Fourier time step: " << vInfo[0].h_gridpoint*vInfo[0].h_gridpoint*.5/nu << ")\n";
	
	const int nsteps = 1;
	
	for(int step=0; step<nsteps; ++step)
	{
		resetOMP(vInfo, *grid);
		//processOMP<Lab,OperatorDiffusionHighOrder>(dt,nu,vInfo,*grid); // this does not work correctly!
		processOMP<Lab,OperatorDiffusion>(dt,nu,vInfo,*grid);
		updateOMP(vInfo, *grid);
		
		time += dt;
	}
	
	stringstream ss;
	ss << path2file << "-bpd" << bpd << ".vti";
	
	dumper.Write(*grid, ss.str());
}

void TestDiffusion::check()
{
	cout << "\tErrors (Linf, L1, L2):\t";
	double Linf = 0.;
	double L1 = 0.;
	double L2 = 0.;
	
	const int size = bpd * FluidBlock::sizeX;
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
#pragma omp parallel for reduction(max:Linf) reduction(+:L1) reduction(+:L2)
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				double p[3];
				info.pos(p, ix, iy);
				
				double error = b(ix, iy).u - _analytical(p[0],p[1],time);
				Linf = max(Linf,abs(error));
				L1 += abs(error);
				L2 += error*error;
			}
	}
	
	L2 = sqrt(L2)/(double)size;
	L1 /= (double)size*size;
	cout << "\t" << Linf << "\t" << L1 << "\t" << L2 << endl;
	myfile << size << " " << Linf << " " << L1 << " " << L2 << endl;
}