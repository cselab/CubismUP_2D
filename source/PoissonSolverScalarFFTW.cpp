//
//  PoissonSolverScalarFFTW.cpp
//  CubismUP_2D
//
//  Created by CSE Lab.
//  Adapted by Christian Conti for non-periodic BC
//  Copyright (c) ETHZ. All rights reserved.
//

#include "PoissonSolverScalarFFTW.h"

int FFTWBase::registered_objects = 0;
bool FFTWBase::initialized = false;

