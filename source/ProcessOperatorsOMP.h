//
//  ProcessOperatorsOMP.h
//  CubismUP_2D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_ProcessOperators_h
#define CubismUP_2D_ProcessOperators_h

#include "Definitions.h"
#include "Shape.h"


// -gradp, divergence, advection
template<typename Lab, typename Kernel>
void processOMP(double dt, vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	BlockInfo * ary = &myInfo.front();
	const int N = myInfo.size();
	
#pragma omp parallel
	{
		Kernel kernel(dt);
		
		Lab mylab;
		mylab.prepare(grid, kernel.stencil_start, kernel.stencil_end, true);
		
#pragma omp for schedule(static)
		for (int i=0; i<N; i++)
		{
			mylab.load(ary[i], 0);
			
			kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
	}
}

// divergence with layer - still useful for diagnostics
template<typename Lab, typename Kernel>
void processOMP(Layer& outputField, vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	BlockInfo * ary = &myInfo.front();
	const int N = myInfo.size();
	
#pragma omp parallel
	{
		Kernel kernel(outputField);
		
		Lab mylab;
		mylab.prepare(grid, kernel.stencil_start, kernel.stencil_end, true);
		
#pragma omp for schedule(static)
		for (int i=0; i<N; i++)
		{
			mylab.load(ary[i], 0, false);
			
			kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
	}
}

// divergence split with layer - still useful for diagnostics
template<typename Lab, typename Kernel>
void processOMP(Layer& outputField, const Real rho0, const Real dt, const int step, vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	BlockInfo * ary = &myInfo.front();
	const int N = myInfo.size();
	
#pragma omp parallel
	{
		Kernel kernel(outputField, rho0, dt, step);
		
		Lab mylab;
		mylab.prepare(grid, kernel.stencil_start, kernel.stencil_end, true);
		
#pragma omp for schedule(static)
		for (int i=0; i<N; i++)
		{
			mylab.load(ary[i], 0, false);
			
			kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
	}
}

double findMaxUOMP(vector<BlockInfo>& myInfo, FluidGrid & grid);
void computeForcesFromVorticity(vector<BlockInfo>& myInfo, FluidGrid & grid, Real ub[2], Real oldAccVort[2], Real rhoS);
#endif
