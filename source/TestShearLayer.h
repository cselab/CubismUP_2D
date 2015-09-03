//
//  TestShearLayer.h
//  CubismUP_2D
//
//  Created by Christian Conti on 4/28/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_2D__TestShearLayer__
#define __CubismUP_2D__TestShearLayer__

#include "Test.h"
#include "GenericCoordinator.h"

class TestShearLayer : public Test
{
protected:
	double nu;
	double dtCFL, dtLCFL, dtFourier;
	double time, endTime;
	const int bpd, bpdRef;
	int step;
	const double rhoS;
	
	// MPI stuff - required for Hypre
	int rank, nprocs;
	
	string path2file;
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
	vector<GenericCoordinator *> pipeline;
	
	FluidGrid * grid;
	FluidGrid * gridRef;
	
	void _ic();
	void _getRefs(const int ix, const int iy, const int ratio, Real &u, Real &v);
	
public:
	TestShearLayer(const int argc, const char ** argv, const int bpd);
	~TestShearLayer();
	
	void run();
	void check();
};

#endif /* defined(__CubismUP_2D__TestShearLayer__) */
