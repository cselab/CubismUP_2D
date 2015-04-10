//
//  CoordinatorCleanTmp.h
//  CubismUP_2D
//
//  Created by Christian Conti on 4/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_CoordinatorCleanTmp_h
#define CubismUP_2D_CoordinatorCleanTmp_h

#include "GenericCoordinator.h"

class CoordinatorCleanTmp : public GenericCoordinator
{
public:
	CoordinatorCleanTmp(FluidGrid * grid) : GenericCoordinator(grid)
	{
	}
	
	void operator()(const double dt)
	{
		check("clean - start");
		/*
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
					b(ix,iy).tmp  = 0;
				}
		}
		*/
		check("clean - end");
	}
	
	string getName()
	{
		return "Clean";
	}
};


#endif
