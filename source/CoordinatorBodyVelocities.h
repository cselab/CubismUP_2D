//
//  CoordinatorBodyVelocities.h
//  CubismUP_2D
//
//  Created by Christian Conti on 3/30/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_CoordinatorComputeBodyVelocities_h
#define CubismUP_2D_CoordinatorComputeBodyVelocities_h

class CoordinatorBodyVelocities : public GenericCoordinator
{
protected:
	Real *uBody, *vBody, *omegaBody;
	Real *lambda;
	Shape * shape;
	
public:
	CoordinatorBodyVelocities(Real * uBody, Real * vBody, Real * omegaBody, Shape * shape, Real * lambda, FluidGrid * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), omegaBody(omegaBody), lambda(lambda), shape(shape)
	{
	}
	
	void operator()(const double dt)
	{
		Real centerOfMass[2];
		shape->getCenterOfMass(centerOfMass);
		
		// gravity acts on the center of gravity (for constant g: center of mass)
		// buoyancy acts on the center of buoyancy (for fully immersed body: center of geometry/centroid)
		double mass = 0;
		double u = 0;
		double v = 0;
		double momOfInertia = 0;
		double angularMomentum = 0;
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static) reduction(+:mass)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			double h = info.h_gridpoint;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				double p[2] = {0,0};
				info.pos(p, ix, iy);
				const double chi = b(ix,iy).chi;
				const double rhochi = b(ix,iy).rho * chi;
				mass += rhochi;
			}
		}
		
#pragma omp parallel for schedule(static) reduction(+:u) reduction(+:v) reduction(+:momOfInertia) reduction(+:angularMomentum)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			double h = info.h_gridpoint;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					double p[2] = {0,0};
					info.pos(p, ix, iy);
					
					double rhochi = b(ix,iy).rho * b(ix,iy).chi;
					u += b(ix,iy).u * rhochi;
					v += b(ix,iy).v * rhochi;
					momOfInertia    += rhochi * ((p[0]-centerOfMass[0])*(p[0]-centerOfMass[0]) + (p[1]-centerOfMass[1])*(p[1]-centerOfMass[1]));
					angularMomentum += rhochi * ((p[0]-centerOfMass[0])*b(ix,iy).v             - (p[1]-centerOfMass[1])*b(ix,iy).u);
				}
		}
	
		*uBody = u / mass;
		*vBody = v / mass;
		*omegaBody = angularMomentum / momOfInertia;
	}
	
	string getName()
	{
		return "BodyVelocities";
	}
};


#endif
