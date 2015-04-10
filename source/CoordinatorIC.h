//
//  CoordinatorIC.h
//  CubismUP_2D
//
//  Created by Christian Conti on 3/27/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_CoordinatorIC_h
#define CubismUP_2D_CoordinatorIC_h

#include "GenericCoordinator.h"
#include "OperatorIC.h"
#include "Shape.h"

class CoordinatorIC : public GenericCoordinator
{
protected:
	Shape * shape;
	const double uinf;
	
public:
	CoordinatorIC(Shape * shape, const double uinf, FluidGrid * grid) : GenericCoordinator(grid), shape(shape), uinf(uinf)
	{
	}
	
	void operator()(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorIC kernel(shape, uinf);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
				kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
		
		check("IC - end");
	}
	
	string getName()
	{
		return "IC";
	}
};

class CoordinatorIC_RT : public GenericCoordinator
{
public:
	CoordinatorIC_RT(FluidGrid * grid) : GenericCoordinator(grid)
	{
	}
	
	void operator()(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorIC_RT kernel;
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
				kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
		
		check("IC - end");
	}
	
	string getName()
	{
		return "IC_RT";
	}
};

#endif
