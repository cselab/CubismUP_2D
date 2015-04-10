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
	double nu;
	double dt;
	
public:
	OperatorDiffusion(double dt, double nu) : nu(nu), dt(dt)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	
	~OperatorDiffusion() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double prefactor = nu * dt / (info.h_gridpoint*info.h_gridpoint);
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				FluidElement& phi = lab(ix,iy);
				FluidElement& phiN = lab(ix,iy+1);
				FluidElement& phiS = lab(ix,iy-1);
				FluidElement& phiE = lab(ix+1,iy);
				FluidElement& phiW = lab(ix-1,iy);
				
				o(ix,iy).tmpU += phi.u + prefactor * (phiN.u + phiS.u + phiE.u + phiW.u - 4.*phi.u);
				o(ix,iy).tmpV += phi.v + prefactor * (phiN.v + phiS.v + phiE.v + phiW.v - 4.*phi.v);
#ifdef _MULTIPHASE_
				o(ix,iy).tmp += phi.rho + prefactor * (phiN.rho + phiS.rho + phiE.rho + phiW.rho - 4.*phi.rho);
#endif // _MULTIPHASE_
			}
	}
};

class OperatorDiffusionHighOrder : public GenericLabOperator
{
private:
	double nu;
	double dt;
	
public:
	OperatorDiffusionHighOrder(double dt, double nu) : nu(nu), dt(dt)
	{
		cout << "This operator needs debugging - it's not working correctly!\n";
		abort();
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = 0;
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 1;
	}
	
	~OperatorDiffusionHighOrder() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double prefactor = nu * dt / (info.h_gridpoint*info.h_gridpoint);
		
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
				
				o(ix,iy).tmpU += phi.u + prefactor/12. * (-(phiN2.u + phiS2.u + phiE2.u + phiW2.u) + 16*(phiN.u + phiS.u + phiE.u + phiW.u) - 60.*phi.u);
				o(ix,iy).tmpV += phi.v + prefactor/12. * (-(phiN2.v + phiS2.v + phiE2.v + phiW2.v) + 16*(phiN.v + phiS.v + phiE.v + phiW.v) - 60.*phi.v);
#ifdef _MULTIPHASE_
				o(ix,iy).tmpV += phi.rho + prefactor/12. * (-(phiN2.rho + phiS2.rho + phiE2.rho + phiW2.rho) + 16*(phiN.rho + phiS.rho + phiE.rho + phiW.rho) - 60.*phi.rho);
#endif // _MULTIPHASE_
			}
	}
};


#endif
