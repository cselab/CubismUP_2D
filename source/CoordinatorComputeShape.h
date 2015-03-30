//
//  CoordinatorComputeShape.h
//  CubismUP_2D
//
//  Created by Christian Conti on 3/30/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_CoordinatorComputeShape_h
#define CubismUP_2D_CoordinatorComputeShape_h

#include "GenericCoordinator.h"
#include "OperatorComputeShape.h"
#include "Shape.h"

class CoordinatorComputeShape : public GenericCoordinator
{
protected:
	Real *uBody, *vBody, *omegaBody;
	Shape * shape;
public:
	CoordinatorComputeShape(Real * uBody, Real * vBody, Real * omegaBody, Shape * shape, FluidGrid * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), omegaBody(omegaBody), shape(shape)
	{
	}
	
	void operator()(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
		Real ub[2] = { *uBody, *vBody };
		shape->updatePosition(ub, *omegaBody, dt);
		
#pragma omp parallel
		{
			OperatorComputeShape kernel(shape);
			
#pragma omp for schedule(static)
			for(int i=0; i<N; i++)
				kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
	}
	
	string getName()
	{
		return "ComputeShape";
	}
};

#endif
