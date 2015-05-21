//
//  InterpolationKernels.h
//  CubismUP_2D
//
//  Created by Christian Conti on 3/6/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_InterpolationKernels_h
#define CubismUP_2D_InterpolationKernels_h

struct Hat
{
	static const int support_start = 0;
	static const int support_end = 2;
	static const int support = 2;
	
	static Real weight(const Real x)
	{
		return 1-abs(x);
	}
};

struct Mp4
{
	static const int support_start = -1;
	static const int support_end = 3;
	static const int support = 4;
	
	static Real weight(const Real x)
	{
		Real absx = abs(x);
		
		// Horner's scheme
		return (absx<2. && absx>=1.) * (((-.5*absx + 2.5)*absx - 4)*absx + 2) + (absx<1.) * (((1.5*absx - 2.5)*absx + 0.)*absx + 1.);
	}
};

struct Ms6
{
	static const int support_start = -2;
	static const int support_end = 4;
	static const int support = 6;
	
	static Real weight(const Real x)
	{
		Real absx = abs(x);
		Real absx2 = absx*absx;
		Real absx3 = absx2*absx;
		Real absx4 = absx2*absx2;
		Real absx5 = absx3*absx2;
		
		if (absx<3 && absx>=2)
			return -1./24.*(absx-2.)*(absx-3.)*(absx-3.)*(absx-3.)*(5*absx-8);
		else if (absx<2 && absx>=1)
			return 1./24.*(absx-1)*(absx-2)*(25*absx3-114*absx2+153*absx-48);
		else if (absx<1)
			return -1./12*(absx-1)*(25*absx4-38*absx3-3*absx2+12*absx+12);
		else
			return 0;
	}
};

#endif
