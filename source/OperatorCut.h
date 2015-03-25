//
//  OperatorCut.h
//  CubismUP_2D
//
//  Created by Christian Conti on 1/10/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_OperatorCut_h
#define CubismUP_2D_OperatorCut_h

class OperatorCut
{
public:
	void operator() (const int width, Layer& v)
	{
		for(int c=0; c<v.nDim; c++)
			for(int iy=0; iy<v.sizeY; iy++)
				for(int ix=0; ix<v.sizeX; ix++)
				{
					if (iy>=width && iy<=v.sizeY-width-1 && ix==width)
					{
						ix= v.sizeX-width-1;
						continue;
					}
					
					v(ix, iy, c) = 0;
				}
	}
};

#endif
