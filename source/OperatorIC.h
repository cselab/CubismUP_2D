//
//  OperatorIC.h
//  CubismUP_2D
//
//	Operates on
//		chi, u, v, rho, p, pOld
//
//  Created by Christian Conti on 1/7/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_OperatorIC_h
#define CubismUP_2D_OperatorIC_h

#include "GenericOperator.h"

class OperatorIC : public GenericOperator
{
private:
	Shape * shape;
	const double uinf;
	
public:
	OperatorIC(Shape * shape, const double uinf) : shape(shape), uinf(uinf) {}
	
    ~OperatorIC() {}
    
    void operator()(const BlockInfo& info, FluidBlock& block) const
    {
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
            for(int ix=0; ix<FluidBlock::sizeX; ++ix)
            {
                Real p[2];
                info.pos(p, ix, iy);
				
				block(ix,iy).u = uinf;
                block(ix,iy).v = 0;
				block(ix,iy).chi = shape->chi(p, info.h_gridpoint);
				
				// assume fluid with density 1
				block(ix,iy).rho = shape->rho(p, info.h_gridpoint);
				
				// this is for testing purposes only! do it the clean way!!
				block(ix,iy).p = 9.81*(1-p[1]);//0;
				block(ix,iy).divU = 9.81*(1-p[1]);//0;
				block(ix,iy).pOld = 0;
            }
    }
};


#endif
