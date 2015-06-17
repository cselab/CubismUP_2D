//
//  CoordinatorAdvection.h
//  CubismUP_2D
//
//  Created by Christian Conti on 3/27/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_CoordinatorAdvection_h
#define CubismUP_2D_CoordinatorAdvection_h

#include "GenericCoordinator.h"
#include "OperatorAdvection.h"
#include <cmath>

template <typename Lab>
class CoordinatorAdvection : public GenericCoordinator
{
protected:
#ifdef _MULTIPHASE_
	Real rhoS;
#endif
	
	inline void reset()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					b(ix,iy).tmpU = 0;
					b(ix,iy).tmpV = 0;
#ifdef _MULTIPHASE_
					b(ix,iy).tmp = 0;
#endif // _MULTIPHASE_
				}
		}
	};
	
	inline void update()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
			for(int i=0; i<N; i++)
			{
				BlockInfo info = vInfo[i];
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;
				
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						b(ix,iy).u = b(ix,iy).tmpU;
						b(ix,iy).v = b(ix,iy).tmpV;
#ifdef _MULTIPHASE_
						//b(ix,iy).chi = b(ix,iy).tmp;
						//b(ix,iy).rho = b(ix,iy).chi * rhoS + (1-b(ix,iy).chi);
						
						// threshold density
#ifdef _PARTICLES_
						Real density = min(max(b(ix,iy).tmp,min((Real)1.,rhoS)),max((Real)1.,rhoS));
						b(ix,iy).rho = density;
#else // _PARTICLES_
	  //b(ix,iy).rho = b(ix,iy).tmp;
						Real density = min(max(b(ix,iy).tmp,min((Real)1.,rhoS)),max((Real)1.,rhoS));
						b(ix,iy).rho = density;
#endif // _PARTICLES_
#endif // _MULTIPHASE_
					}
			}
	}
	
	inline void advect(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
#ifndef _PARTICLES_
			// this is wrong - using -u instead of u?
			OperatorAdvectionUpwind3rdOrder kernel(dt);
			//OperatorAdvectionFD kernel(dt);
			
			Lab mylab;
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
#else // _PARTICLES_
			//OperatorAdvection<Hat> kernel(dt);
			//OperatorAdvection<Lambda2> kernel(dt);
			//OperatorAdvection<Mp4> kernel(dt);
			OperatorAdvection<Ms6> kernel(dt);
			
			Lab mylab;
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
#endif // _PARTICLES_
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
public:
#ifndef _MULTIPHASE_
	CoordinatorAdvection(FluidGrid * grid) : GenericCoordinator(grid)
#else
	CoordinatorAdvection(FluidGrid * grid, Real rhoS) : GenericCoordinator(grid), rhoS(rhoS)
#endif
	{
	}
	
	void operator()(const double dt)
	{
		check("advection - start");

		reset();
		advect(dt);
		update();
		
		check("advection - end");
	}
	
	string getName()
	{
		return "Advection";
	}
};

template <typename Lab>
class CoordinatorTransport : public GenericCoordinator
{
protected:
	inline void reset()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					b(ix,iy).tmp = 0;
		}
	};
	
	inline void update()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					b(ix,iy).rho = b(ix,iy).tmp;
		}
	}
	
public:
	CoordinatorTransport(FluidGrid * grid) : GenericCoordinator(grid)
	{
	}
	
	void operator()(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
		reset();
		
#pragma omp parallel
		{
			//OperatorTransport<Hat> kernel(dt);
			//OperatorTransport<Mp4> kernel(dt);
			OperatorTransport<Ms6> kernel(dt);
			
			Lab mylab;
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
		
		update();
	}
	
	string getName()
	{
		return "Transport";
	}
};

#endif
