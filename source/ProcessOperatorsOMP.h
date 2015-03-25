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
		for(int i=0; i<N; i++)
		{
			mylab.load(ary[i], 0);
			
			kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
	}
}

// diffusion
template<typename Lab, typename Kernel>
void processOMP(double dt, double coeff, vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	BlockInfo * ary = &myInfo.front();
	const int N = myInfo.size();
	
#pragma omp parallel
	{
		Kernel kernel(dt, coeff);
		
		Lab mylab;
		mylab.prepare(grid, kernel.stencil_start, kernel.stencil_end, false);
		
#pragma omp for schedule(static)
		for(int i=0; i<N; i++)
		{
			mylab.load(ary[i], 0);
			kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
	}
}

// split
template<typename Lab, typename Kernel>
void processOMP(double dt, double rho0, int step, vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	BlockInfo * ary = &myInfo.front();
	const int N = myInfo.size();
	
#pragma omp parallel
	{
		Kernel kernel(dt, rho0, step);
		
		Lab mylab;
		mylab.prepare(grid, kernel.stencil_start, kernel.stencil_end, false);
		
#pragma omp for schedule(static)
		for(int i=0; i<N; i++)
		{
			mylab.load(ary[i], 0);
			kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
	}
}

// penalization
template<typename Kernel>
void processOMP(double dt, Real uBody, Real vBody, Real omegaBody, Real xCenterOfMass, Real yCenterOfMass, double lambda, vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	BlockInfo * ary = &myInfo.front();
	const int N = myInfo.size();
	
#pragma omp parallel
	{
		Kernel kernel(dt, uBody, vBody, omegaBody, xCenterOfMass, yCenterOfMass, lambda);
		
#pragma omp for schedule(static)
		for(int i=0; i<N; i++)
			kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
	}
}

// Chi field computation
template<typename Kernel>
void processOMP(Shape * shape, vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	BlockInfo * ary = &myInfo.front();
	const int N = myInfo.size();
	
#pragma omp parallel
	{
		Kernel kernel(shape);
		
#pragma omp for schedule(static)
		for(int i=0; i<N; i++)
			kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
	}
}

// gravity field computation - with hydrostatic pressure term
template<typename Kernel>
void processOMP_hydrostaticTerm(Real v[2], double factor, double rhoS, vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	BlockInfo * ary = &myInfo.front();
	const int N = myInfo.size();
	
	Real volS = 0;
	Real volF = 0;
	
#pragma omp parallel for schedule(static) reduction(+:volS) reduction(+:volF)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		Real dh = info.h_gridpoint;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				volS += b(ix,iy).chi*dh*dh;
				volF += (1-b(ix,iy).chi)*dh*dh;
			}
	}
	
	Real hydrostaticFactor = (rhoS*volS+1.*volF)/(volS+volF);
	
#pragma omp parallel
	{
		Kernel kernel(v, factor, hydrostaticFactor);
		
#pragma omp for schedule(static)
		for(int i=0; i<N; i++)
			kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
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
		for(int i=0; i<N; i++)
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
		for(int i=0; i<N; i++)
		{
			mylab.load(ary[i], 0, false);
			
			kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
	}
}

void updateOMP(vector<BlockInfo>& myInfo, FluidGrid & grid);
void updatePressuresOMP(vector<BlockInfo>& myInfo, FluidGrid & grid);
void resetOMP(vector<BlockInfo>& myInfo, FluidGrid & grid);
double findMaxUOMP(vector<BlockInfo>& myInfo, FluidGrid & grid);
void computeBodyVelocity(vector<BlockInfo>& myInfo, FluidGrid & grid, Real ub[2], Real& angularU, Real rhoS, Real g[2], Real dt, Real lambda);
void computeForcesFromVorticity(vector<BlockInfo>& myInfo, FluidGrid & grid, Real ub[2], Real oldAccVort[2], Real rhoS);
#endif
