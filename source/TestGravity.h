//
//  TestGravity.h
//  CubismUP_2D
//
//  Created by Christian Conti on 2/3/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_2D__TestGravity__
#define __CubismUP_2D__TestGravity__

#include <stdio.h>
#include "Test.h"
#include "OperatorGravity.h"

class TestGravity : public Test
{
private:
	double time;
	int bpd;
	
	Real gravity[2];
	
	string path2file;
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
	FluidGrid * grid;
	
	void _ic();
	
public:
	TestGravity(const int argc, const char ** argv, const int bpd);
	~TestGravity();
	
	void run();
	void check();
};

#endif /* defined(__CubismUP_2D__TestGravity__) */
