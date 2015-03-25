//
//  TestTranslation.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 3/16/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestTranslation.h"
#include "ProcessOperatorsOMP.h"
#include "OperatorComputeShape.h"
#include "OperatorVorticity.h"
#include <sstream>

void TestTranslation::_ic()
{
	Real center[2] = {.25,.25};
	Real semiAxis[2] = {.05,.15};
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	bool bPeriodic[2] = {false,false};
	shape = new Ellipse(center, semiAxis, (Real)M_PI/4, (Real)2, (Real)2, (Real)2, bPeriodic);
	
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
		for(int ix=0; ix<FluidBlock::sizeX; ix++)
		{
			Real p[2];
			info.pos(p, ix, iy);
			
			// this is for testCase==1, no effect on testCase==0
			b(ix,iy).u = 1;
			b(ix,iy).v = 1;
			
			b(ix,iy).chi = shape->chi(p, info.h_gridpoint);
			
			// assume fluid with density 1
			b(ix,iy).rho = shape->rho(p, info.h_gridpoint);
			
			// this is for testing purposes only! do it the clean way!!
			b(ix,iy).p = 0;
			b(ix,iy).divU = 0;
			b(ix,iy).pOld = 0;
		}
	}
}

TestTranslation::TestTranslation(const int argc, const char ** argv, const int testCase, const int bpd) : Test(argc,argv), testCase(testCase), bpd(bpd)
{
	grid = new FluidGrid(bpd,bpd,1);
	
	path2file = parser("-file").asString("../data/testTranslation");
	_ic();
	
}

TestTranslation::~TestTranslation()
{
	delete grid;
}

void TestTranslation::run()
{
	const int sizeX = bpd * FluidBlock::sizeX;
	const int sizeY = bpd * FluidBlock::sizeY;
	
	for (int step=0; step<50; step++)
	{
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		
		Real u[2] = {0,0};
		Real omega = 0;
		const Real dt = 1e-2;
		
		if (testCase==0)
		{
			u[0] = 1;
			u[1] = 1;
		}
		else if (testCase==1)
		{
			Real rhoS = shape->getRhoS();
			Real gravity[2] = {0,0};
			Real lambda = 1;
			computeBodyVelocity(vInfo, *grid, u, omega, rhoS, gravity, dt, lambda);
		}
		else
			abort();
		
		shape->updatePosition(u,omega,dt);
		processOMP<OperatorComputeShape>(shape, vInfo, *grid);
		
		// dump
		if (step%10==0)
		{
			stringstream ss;
			ss << path2file << "-" << bpd << "-" << step << ".vti";
			
			dumper.Write(*grid, ss.str());
			
			Layer vorticity(sizeX,sizeY,1);
			processOMP<Lab, OperatorVorticity>(vorticity,vInfo,*grid);
			stringstream sVort;
			sVort << path2file << "Vorticity-" << bpd << "-" << step << ".vti";
			dumpLayer2VTK(step,sVort.str(),vorticity,1);
		}
	}
}

void TestTranslation::check()
{
	// the only thing to check here is the orientation
	Real p[2];
	shape->getPosition(p);
	cout << "Translation error X: " << p[0] - .75 << endl;
	cout << "Translation error Y: " << p[1] - .75 << endl;
}