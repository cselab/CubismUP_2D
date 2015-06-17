//
//  OperatorDiffusion.h
//  CubismUP_2D
//
//	Operates on
//		tmpU, tmpV
//
//  Created by Christian Conti on 1/7/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_OperatorDiffusion_h
#define CubismUP_2D_OperatorDiffusion_h

#include "GenericOperator.h"

class OperatorDiffusion : public GenericLabOperator
{
private:
	const double mu;
	double dt;
	const int stage;
	
public:
	OperatorDiffusion(double dt, double mu, const int stage) : mu(mu), dt(dt), stage(stage)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
		
		//dt = (stage==0) ? dt*.5 : dt;
	}
	
	~OperatorDiffusion() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double prefactor = mu * dt / (info.h_gridpoint*info.h_gridpoint);
		
		// stage 1 of RK2
		if (stage==0)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					FluidElement& phi = lab(ix,iy);
					FluidElement& phiN = lab(ix,iy+1);
					FluidElement& phiS = lab(ix,iy-1);
					FluidElement& phiE = lab(ix+1,iy);
					FluidElement& phiW = lab(ix-1,iy);
					
					o(ix,iy).tmpU = phi.u + prefactor/phi.rho * (phiN.u + phiS.u + phiE.u + phiW.u - phi.u*4.);
					o(ix,iy).tmpV = phi.v + prefactor/phi.rho * (phiN.v + phiS.v + phiE.v + phiW.v - phi.v*4.);
					//o(ix,iy).tmpU = phi.u + prefactor * (phiN.u + phiS.u + phiE.u + phiW.u - phi.u*4.);
					//o(ix,iy).tmpV = phi.v + prefactor * (phiN.v + phiS.v + phiE.v + phiW.v - phi.v*4.);
				}
		/*
		// stage 2 of RK2
		else if (stage==1)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					FluidElement& phi = lab(ix,iy);
					FluidElement& phiN = lab(ix,iy+1);
					FluidElement& phiS = lab(ix,iy-1);
					FluidElement& phiE = lab(ix+1,iy);
					FluidElement& phiW = lab(ix-1,iy);
					
					o(ix,iy).tmpU = phi.u + prefactor * (phiN.tmpU + phiS.tmpU + phiE.tmpU + phiW.tmpU - 4.*phi.tmpU);
					o(ix,iy).tmpV = phi.v + prefactor * (phiN.tmpV + phiS.tmpV + phiE.tmpV + phiW.tmpV - 4.*phi.tmpV);
				}
		 */
	}
};

class OperatorDiffusionHighOrder : public GenericLabOperator
{
private:
	double mu;
	double dt;
	const int stage;
	
public:
	OperatorDiffusionHighOrder(double dt, double mu, const int stage) : mu(mu), dt(dt), stage(stage)
	{
		// it might be working - but the error is too small to measure!
		//cout << "This operator needs debugging - it's not working correctly!\n";
		//abort();
		
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = 0;
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 1;
		
		//dt = (stage==0) ? dt*.5 : dt;
	}
	
	~OperatorDiffusionHighOrder() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double prefactor = mu * dt / (12*info.h_gridpoint*info.h_gridpoint);
		
		if (stage==0)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					FluidElement& phi = lab(ix,iy);
					FluidElement& phiN = lab(ix,iy+1);
					FluidElement& phiS = lab(ix,iy-1);
					FluidElement& phiE = lab(ix+1,iy);
					FluidElement& phiW = lab(ix-1,iy);
					FluidElement& phiN2 = lab(ix,iy+2);
					FluidElement& phiS2 = lab(ix,iy-2);
					FluidElement& phiE2 = lab(ix+2,iy);
					FluidElement& phiW2 = lab(ix-2,iy);
				
					o(ix,iy).tmpU += phi.u + prefactor * (-(phiN2.u + phiS2.u + phiE2.u + phiW2.u) + 16*(phiN.u + phiS.u + phiE.u + phiW.u) - 60.*phi.u);
					o(ix,iy).tmpV += phi.v + prefactor * (-(phiN2.v + phiS2.v + phiE2.v + phiW2.v) + 16*(phiN.v + phiS.v + phiE.v + phiW.v) - 60.*phi.v);
				}
	}
};


#endif
