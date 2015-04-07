//
//  CoordinatorPenalization.h
//  CubismUP_2D
//
//  Created by Christian Conti on 3/27/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_CoordinatorPenalization_h
#define CubismUP_2D_CoordinatorPenalization_h

#include "GenericCoordinator.h"
#include "OperatorPenalization.h"
#include "Shape.h"

class CoordinatorPenalizationFixed : public GenericCoordinator
{
protected:
	Shape * shape;
	Real lambda;
	
public:
	CoordinatorPenalizationFixed(Shape * shape, Real lambda, FluidGrid * grid) : GenericCoordinator(grid), shape(shape), lambda(lambda)
	{
	}
	
	void operator()(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Real com[2];
			shape->getPosition(com);
			OperatorPenalization kernel(dt, 0, 0, 0, com[0], com[1], lambda);
			
#pragma omp for schedule(static)
			for(int i=0; i<N; i++)
				kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
	}
	
	string getName()
	{
		return "Penalization";
	}
};

class CoordinatorPenalization : public GenericCoordinator
{
protected:
	Real *uBody, *vBody, *omegaBody;
	Shape * shape;
	Real lambda;
	
public:
	CoordinatorPenalization(Real * uBody, Real * vBody, Real * omegaBody, Shape * shape, Real lambda, FluidGrid * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), omegaBody(omegaBody), shape(shape), lambda(lambda)
	{
	}
	
	void operator()(const double dt)
	{
		check("penalization - start");
		
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Real com[2];
			shape->getPosition(com);
			OperatorPenalization kernel(dt, *uBody, *vBody, *omegaBody, com[0], com[1], lambda);
			
#pragma omp for schedule(static)
			for(int i=0; i<N; i++)
				kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
		
		check("penalization - end");
	}
	
	string getName()
	{
		return "Penalization";
	}
};

#endif
