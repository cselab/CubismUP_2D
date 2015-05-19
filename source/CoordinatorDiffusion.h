//
//  CoordinatorDiffusion.h
//  CubismUP_2D
//
//  Created by Christian Conti on 3/27/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_CoordinatorDiffusion_h
#define CubismUP_2D_CoordinatorDiffusion_h

#include "GenericCoordinator.h"
#include "OperatorDiffusion.h"

template <typename Lab>
class CoordinatorDiffusion : public GenericCoordinator
{
protected:
	const double coeff;
	
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
					/*
					 // there is no diffusion of density!
#ifdef _MULTIPHASE_
					b(ix,iy).tmp = 0;
#endif // _MULTIPHASE_
					 */
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
					/*
					 // there is no diffusion of density!
#ifdef _MULTIPHASE_
					b(ix,iy).tmp = -b(ix,iy).tmp/2;
#endif // _MULTIPHASE_
					 */
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
					/*
					 // there is no diffusion of density!
#ifdef _MULTIPHASE_
					b(ix,iy).rho = b(ix,iy).tmp;
#endif // _MULTIPHASE_
					 */
				}
		}
	 }
	
	inline void diffuse(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorDiffusion kernel(dt, coeff);
			//OperatorDiffusionHighOrder  kernel(dt, coeff);
			
			Lab mylab;
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
public:
	CoordinatorDiffusion(const double coeff, FluidGrid * grid) : GenericCoordinator(grid), coeff(coeff)
	{
	}
	
	void operator()(const double dt)
	{
		check("diffusion - start");
		
		reset();
		diffuse(dt);
		update();
		/*
		resetHeun();
		diffuse(dt/2);
		update();
		*/
		check("diffusion - end");
	}
	
	string getName()
	{
		return "Diffusion";
	}
};
#endif
