//
//  TestPenalization.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 2/3/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestPenalization.h"
#include <sstream>
#include <cmath>

#include "CoordinatorIC.h"
#include "CoordinatorPenalization.h"

void TestPenalization::_ic()
{
	Real center[2] = {.5,.5};
	Real radius = .1;
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	bool bPeriodic[2] = {false,false};
	shape = new Disk(center, radius, (Real).1, (Real)2, (Real)2, bPeriodic);
	
	CoordinatorIC coordIC(shape,1.,grid);
	coordIC(0);
}

TestPenalization::TestPenalization(const int argc, const char ** argv, const int bpd) : Test(argc, argv), bpd(bpd), lambda(1e4), uBody{0,0.5}
{
	grid = new FluidGrid(bpd,bpd,1);
	
	// setup initial condition
	
	// output settings
	path2file = parser("-file").asString("../data/testPenalization");
	_ic();
}

TestPenalization::~TestPenalization()
{
	delete grid;
	delete shape;
}

void TestPenalization::run()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	const double dt = 1e-4;
	Real omegaBody = 0;
	CoordinatorPenalization coordPenalization(&uBody[0], &uBody[1], &omegaBody, shape, lambda, grid);
	
	for (int i=0; i<10; i++)
	{
		//Real centerOfMass[2] = {0,0};
		//shape->getPosition(centerOfMass);
		//processOMP<OperatorPenalization>(dt,uBody[0],uBody[1],omegaBody,centerOfMass[0],centerOfMass[1],lambda,vInfo,*grid);
		coordPenalization(dt);
	
		stringstream ss;
		ss << path2file << "-bpd" << bpd << "-step" << i << ".vti";
	
		dumper.Write(*grid, ss.str());
	}
}

void TestPenalization::check()
{
}
