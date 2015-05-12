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

template <typename Lab>
class CoordinatorAdvection : public GenericCoordinator
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
				{
					b(ix,iy).tmpU = 0;
					b(ix,iy).tmpV = 0;
#ifdef _MULTIPHASE_
					b(ix,iy).tmp = 0;
#endif // _MULTIPHASE_
				}
		}
	};
	
	inline void resetHeun()
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
					b(ix,iy).tmpU = -b(ix,iy).tmpU/2.;
					b(ix,iy).tmpV = -b(ix,iy).tmpV/2.;
#ifdef _MULTIPHASE_
					b(ix,iy).tmp = -b(ix,iy).tmp/2;
#endif // _MULTIPHASE_
				}
		}
	};
	
	inline void update(const int stage=1)
	{
		const int N = vInfo.size();
		
		if (stage==0)
		{
#pragma omp parallel for schedule(static)
			for(int i=0; i<N; i++)
			{
				BlockInfo info = vInfo[i];
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;
				
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						b(ix,iy).rk2u = b(ix,iy).tmpU;
						b(ix,iy).rk2v = b(ix,iy).tmpV;
						
						// no need to do anything with density
					}
			}
		}
		else if (stage==1)
		{
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
						b(ix,iy).rho = b(ix,iy).tmp;
#endif // _MULTIPHASE_
					}
			}
		}
	}
	
	inline void advect(const double dt, const int stage)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
#ifndef _PARTICLES_
			// this is wrong - using -u instead of u?
			OperatorAdvectionFD kernel(dt, stage);
#else // _PARTICLES_
			OperatorAdvection<Mp4> kernel(dt, stage);
			//OperatorAdvection<Ms6> kernel(dt, stage);
#endif // _PARTICLES_
			
			Lab mylab;
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
public:
	CoordinatorAdvection(FluidGrid * grid) : GenericCoordinator(grid)
	{
	}
	
	void operator()(const double dt)
	{
		check("advection - start");
		
		// Euler
		//*
		reset();
		advect(dt,0);
		update(1);
		//*/
		
		// midpoint
		/*
		reset();
		advect(dt/2,0);
		update(0);
		reset();
		advect(dt,1);
		update(1);
		//*/
		
		/*
		 // check!
		reset();
		advect(dt);
		update();
		// requires p2m
		resetHeun();
		advect(dt/2);
		update();
		//*/
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
	/*
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
	 */
	
public:
	CoordinatorTransport(FluidGrid * grid) : GenericCoordinator(grid)
	{
	}
	
	void operator()(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
		//reset();
		
#pragma omp parallel
		{
			OperatorTransport<Mp4> kernel(dt);
			
			Lab mylab;
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
		
		//update();
	}
	
	string getName()
	{
		return "Transport";
	}
};

#endif
