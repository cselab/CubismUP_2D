//
//  Simulation_Fluid.h
//  CubismUP_2D
//
//	Base class for fluid simulations from which any fluid simulation case should inherit
//	Contains the base structure and interface that any fluid simulation class should have
//
//  Created by Christian Conti on 3/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_Simulation_Fluid_h
#define CubismUP_2D_Simulation_Fluid_h

#include "Definitions.h"
#include "ProcessOperatorsOMP.h"
#include "OperatorVorticity.h"
#include "GenericCoordinator.h"
#include "GenericOperator.h"

#include <vector>

class Simulation_Fluid
{
protected:
	ArgumentParser parser;
	
	vector<GenericCoordinator *> pipeline;
	
	// grid
	int bpdx, bpdy;
	FluidGrid * grid;
	
	// simulation status
	int step, nsteps;
	double dt, time, endTime;
	
	// verbose
	bool verbose;
	
	// output
	int dumpFreq;
	double dumpTime;
	string path2file;
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
	virtual void _diagnostics() = 0;
	virtual void _ic() = 0;
	virtual double _nonDimensionalTime() = 0;
	
	virtual void _dump(double & nextDumpTime)
	{
		const int sizeX = bpdx * FluidBlock::sizeX;
		const int sizeY = bpdy * FluidBlock::sizeY;
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		
		if((dumpFreq>0 && step % dumpFreq == 0) || (dumpTime>0 && abs(_nonDimensionalTime()-nextDumpTime) < 10*std::numeric_limits<Real>::epsilon()))
		{
			nextDumpTime += dumpTime;
			
			stringstream ss;
			ss << path2file << "-" << step << ".vti";
			cout << ss.str() << endl;
			
			dumper.Write(*grid, ss.str());
			
			Layer vorticity(sizeX,sizeY,1);
			processOMP<Lab, OperatorVorticity>(vorticity,vInfo,*grid);
			stringstream sVort;
			sVort << path2file << "Vorticity-" << step << ".vti";
			dumpLayer2VTK(step,sVort.str(),vorticity,1);
		}
	}
	
public:
	Simulation_Fluid(const int argc, const char ** argv) : parser(argc,argv), step(0), time(0), dt(0)
	{
		// initialize grid
		parser.set_strict_mode();
		bpdx = parser("-bpdx").asInt();
		bpdy = parser("-bpdy").asInt();
		grid = new FluidGrid(bpdx,bpdy,1);
		assert(grid != NULL);
		
		// simulation ending parameters
		parser.unset_strict_mode();
		nsteps = parser("-nsteps").asInt(0);		// nsteps==0   means that this stopping criteria is not active
		endTime = parser("-tend").asDouble(0);		// endTime==0  means that this stopping criteria is not active
		
		// output parameters
		dumpFreq = parser("-fdump").asDouble(0);	// dumpFreq==0 means that this dumping frequency (in #steps) is not active
		dumpTime = parser("-tdump").asDouble(0);	// dumpTime==0 means that this dumping frequency (in time)   is not active
		path2file = parser("-file").asString("../data/Simulation_Fluid");
		
		verbose = parser("-verbose").asBool(false);
	}
	
	virtual ~Simulation_Fluid()
	{
		delete grid;
		
		while(!pipeline.empty())
		{
			GenericCoordinator * g = pipeline.back();
			pipeline.pop_back();
			delete g;
		}
	}
	
	virtual void simulate() = 0;
};

#endif
