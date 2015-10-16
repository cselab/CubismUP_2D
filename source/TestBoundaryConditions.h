//
//  TestBoundaryConditions.h
//  CubismUP_2D
//
//  Created by Christian Conti on 10/16/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_2D__TestBoundaryConditions__
#define __CubismUP_2D__TestBoundaryConditions__

#include <stdio.h>
#include "Test.h"

class TestBoundaryConditions : public Test
{
private:
    int offset;
    
    string path2file;
    SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
    
    FluidGrid * grid;
    
    void _ic();
    
public:
    TestBoundaryConditions(const int argc, const char ** argv);
    ~TestBoundaryConditions();
    
    void run();
    void check();
};

#endif /* defined(__CubismUP_2D__TestBoundaryConditions__) */
