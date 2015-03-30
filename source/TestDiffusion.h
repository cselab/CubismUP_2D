//
//  TestDiffusion.h
//  CubismUP_2D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_2D__TestDiffusion__
#define __CubismUP_2D__TestDiffusion__

#include <stdio.h>
#include "Test.h"
//#include "OperatorDiffusion.h"

class TestDiffusion : public Test
{
private:
	double nu;
    double time;
	const int bpd;
    
    string path2file;
    SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
    
    FluidGrid * grid;
    
    void _ic();
    double _analytical(double ix, double iy, double t);
    
public:
    TestDiffusion(const int argc, const char ** argv, const int bpd);
	~TestDiffusion();
		
    void run();
    void check();
};

#endif /* defined(__CubismUP_2D__TestDiffusion__) */
