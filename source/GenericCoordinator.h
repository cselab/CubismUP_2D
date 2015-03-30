//
//  GenericCoordinator.h
//  CubismUP_2D
//
//	This class serves as the interface for a coordinator object
//	A coordinator object schedules the processing of blocks with its operator
//
//  Created by Christian Conti on 3/27/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_GenericCoordinator_h
#define CubismUP_2D_GenericCoordinator_h

#include "Definitions.h"

class GenericCoordinator
{
protected:
	FluidGrid * grid;
	vector<BlockInfo> vInfo;
	
public:
	GenericCoordinator(FluidGrid * grid) : grid(grid)
	{
		vInfo = grid->getBlocksInfo();
	}
	
	virtual void operator()(const double dt) = 0;
	
	virtual string getName() = 0;
};

#endif
