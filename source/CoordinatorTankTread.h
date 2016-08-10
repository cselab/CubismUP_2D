//
//  CoordinatorTankTread.h
//  CubismUP_2D
//
//  Created by Christian Conti on 12/14/15.
//  Copyright © 2015 ETHZ. All rights reserved.
//

#ifndef CoordinatorTankTread_h
#define CoordinatorTankTread_h

class CoordinatorTankTread : public GenericCoordinator
{
protected:
	Shape * shape;
	
public:
	CoordinatorTankTread(Shape * shape, FluidGrid * grid) : GenericCoordinator(grid), shape(shape)
	{
	}
	
	void operator()(const double dt)
	{
#ifdef _TANKTREADING_
		Real centerOfMass[2];
		shape->getCenterOfMass(centerOfMass);
		Real omega = shape->getSteering();
		
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			double h = info.h_gridpoint;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					double p[2];
					info.pos(p, ix, iy);
					p[0] -= centerOfMass[0];
					p[1] -= centerOfMass[1];
					
					const Real band = 4. * b(ix,iy).chi * (1.-b(ix,iy).chi); // the factor 4 is needed to make the peak=1 and it comes from the vertical shift of the cosines
					// alternative approach: use interpolation schemes from Immersed Boundaries (?)
					b(ix,iy).u -= band * p[1] * omega;
					b(ix,iy).v += band * p[0] * omega;
				}
		}
#endif
	}
	
	string getName()
	{
		return "TankTread";
	}
};

#endif /* CoordinatorTankTread_h */
