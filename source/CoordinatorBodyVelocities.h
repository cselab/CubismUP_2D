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
	Real lambda;
	
public:
	CoordinatorBodyVelocities(Real * uBody, Real * vBody, Real * omegaBody, const Real lambda, FluidGrid * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), omegaBody(omegaBody), lambda(lambda)
	{
	}
	
	void operator()(const double dt)
	{
		double centerTmpX = 0;
		double centerTmpY = 0;
		double mass = 0;
		double u = 0;
		double v = 0;
		double momOfInertia = 0;
		double angularMomentum = 0;
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static) reduction(+:centerTmpX) reduction(+:centerTmpY) reduction(+:mass)
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
					centerTmpX += p[0] * rhochi;
					centerTmpY += p[1] * rhochi;
					mass += rhochi;
				}
		}
		
		centerTmpX /= mass;
		centerTmpY /= mass;
		
		//*
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
					momOfInertia    += rhochi * ((p[0]-centerTmpX)*(p[0]-centerTmpX) + (p[1]-centerTmpY)*(p[1]-centerTmpY));
					angularMomentum += rhochi * ((p[0]-centerTmpX)*b(ix,iy).v        - (p[1]-centerTmpY)*b(ix,iy).u);
				}
		}
		
		*uBody = u / mass;
		*vBody = v / mass;
		*omegaBody = angularMomentum / momOfInertia;
		
		/*/
		 #pragma omp parallel for schedule(static) reduction(+:u) reduction(+:v)
		 for(int i=0; i<N; i++)
		 {
		 BlockInfo info = vInfo[i];
		 FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		 
		 Real h = info.h_gridpoint;
		 
		 for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
		 u += (b(ix,iy).u-ub[0]) * b(ix,iy).chi;
		 v += (b(ix,iy).v-ub[1]) * b(ix,iy).chi;
			}
		 }
		 
		 ub[0] += dt*u*lambda / mass;
		 ub[1] += dt*v*lambda / mass;
		 //*/
	}
	
	string getName()
	{
		return "BodyVelocities";
	}
};


#endif