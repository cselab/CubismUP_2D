//
//  common.h
//  CubismUP_2D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_common_h
#define CubismUP_2D_common_h

#include <cassert>
#include <sstream>
#include <cmath>
#include <cstdio>

// utmost import to be defined before including cubism

// double precision fftw not installed currently
#ifndef _SP_COMP_
typedef double Real;
#else // _SP_COMP_
typedef float Real;
#endif // _SP_COMP_

//this is all cubism file we need
#include <ArgumentParser.h>
#include <Grid.h>
#include <BlockInfo.h>
#include <SerializerIO_ImageVTK.h>
#include <BlockProcessing.h>
#include <BlockLab.h>
#ifdef _MULTIGRID_
#include <mpi.h>
#endif // _MULTIGRID_

#include "Timer.h"

#endif
