//
//  OperatorAdvection.h
//  CubismUP_2D
//
//	Operates on
//		tmpU, tmpV
//
//  Created by Christian Conti on 1/7/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_OperatorAdvection_h
#define CubismUP_2D_OperatorAdvection_h

#include <cmath>
#include "InterpolationKernels.h"
#include "GenericOperator.h"

template <typename RemeshingKernel>
class OperatorAdvection : public GenericLabOperator
{
private:
	double dt;
	
public:
	OperatorAdvection(double dt) : dt(dt)
	{
		stencil_start[0] = RemeshingKernel::support_start-2;
		stencil_start[1] = RemeshingKernel::support_start-2;
		stencil_start[2] = 0;
		
		stencil_end[0] = RemeshingKernel::support_end+1;
		stencil_end[1] = RemeshingKernel::support_end+1;
		stencil_end[2] = 1;
	}
	~OperatorAdvection() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double dh = info.h_gridpoint;
		const double invdh = 1./dh;
		
		const int support_start = RemeshingKernel::support_start;
		const int support_end   = RemeshingKernel::support_end;
		
		const int bx = info.index[0]*FluidBlock::sizeX;
		const int by = info.index[1]*FluidBlock::sizeY;
		
		for(int iy=stencil_start[1]; iy<FluidBlock::sizeY+stencil_end[1]-1; ++iy)
			for(int ix=stencil_start[0]; ix<FluidBlock::sizeX+stencil_end[0]-1; ++ix)
			{
				double p[2];
				info.pos(p,ix,iy);
				
				FluidElement& particle = lab(ix,iy);
				
				p[0] += dt * particle.u;
				p[1] += dt * particle.v;
				
				// nearest point with lower index
#ifndef _VERTEXCENTERED_
				const double px = p[0]*invdh-.5;
				const double py = p[1]*invdh-.5;
#else
				const double px = p[0]*invdh;
				const double py = p[1]*invdh;
#endif
				const int fpx = (int)floor(px);
				const int fpy = (int)floor(py);
				
				// compute weights
				double wx[RemeshingKernel::support], wy[RemeshingKernel::support];
				for (int i=support_start; i<support_end; i++)
				{
					wx[i-support_start] = RemeshingKernel::weight(px-(fpx+i));
					wy[i-support_start] = RemeshingKernel::weight(py-(fpy+i));
				}
				
				// scatter only to elements within the block, elements outside the block are taken care by other blocks
				for (int j=support_start; j<support_end; j++)
					for (int i=support_start; i<support_end; i++)
					{
						if (fpx+i>=bx && fpx+i<bx+FluidBlock::sizeX &&
							fpy+j>=by && fpy+j<by+FluidBlock::sizeY)
						{
							const int lfpx = fpx+i - bx;
							const int lfpy = fpy+j - by;
							const double weight = wx[i-support_start] * wy[j-support_start];
							o(lfpx,lfpy).tmpU   += weight * particle.u;
							o(lfpx,lfpy).tmpV   += weight * particle.v;
#ifdef _MULTIPHASE_
							o(lfpx,lfpy).tmp += weight * particle.rho;
#endif // _MULTIPHASE_
						}
					}
			}
	}
};

template <typename RemeshingKernel>
class OperatorTransport : public GenericLabOperator
{
private:
	double dt;
	
public:
	OperatorTransport(double dt) : dt(dt)
	{
		stencil_start[0] = RemeshingKernel::support_start-2;
		stencil_start[1] = RemeshingKernel::support_start-2;
		stencil_start[2] = 0;
		
		stencil_end[0] = RemeshingKernel::support_end+1;
		stencil_end[1] = RemeshingKernel::support_end+1;
		stencil_end[2] = 1;
	}
	~OperatorTransport() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double dh = info.h_gridpoint;
		const double invdh = 1./dh;
		
		// particles
		const int support_start = RemeshingKernel::support_start;
		const int support_end   = RemeshingKernel::support_end;
		
		const int bx = info.index[0]*FluidBlock::sizeX;
		const int by = info.index[1]*FluidBlock::sizeY;
		
		for(int iy=stencil_start[1]; iy<FluidBlock::sizeY+stencil_end[1]-1; ++iy)
			for(int ix=stencil_start[0]; ix<FluidBlock::sizeX+stencil_end[0]-1; ++ix)
			{
				double p[2];
				info.pos(p,ix,iy);
				
				FluidElement& particle = lab(ix,iy);
				
				p[0] += dt * particle.u;
				p[1] += dt * particle.v;
				
				// nearest point with lower index
				const double px = p[0]*invdh-.5;
				const double py = p[1]*invdh-.5;
				const int fpx = (int)floor(px);
				const int fpy = (int)floor(py);
				
				// compute weights
				double wx[RemeshingKernel::support], wy[RemeshingKernel::support];
				for (int i=support_start; i<support_end; i++)
				{
					wx[i-support_start] = RemeshingKernel::weight(px-(fpx+i));
					wy[i-support_start] = RemeshingKernel::weight(py-(fpy+i));
				}
				
				// scatter only to elements within the block, elements outside the block are taken care by other blocks
				for (int j=support_start; j<support_end; j++)
					for (int i=support_start; i<support_end; i++)
					{
						if (fpx+i>=bx && fpx+i<bx+FluidBlock::sizeX &&
							fpy+j>=by && fpy+j<by+FluidBlock::sizeY)
						{
							const int lfpx = fpx+i - bx;
							const int lfpy = fpy+j - by;
							const double weight = wx[i-support_start] * wy[j-support_start];
							o(lfpx,lfpy).tmp += weight * particle.rho;
						}
					}
			}
	}
};

class OperatorAdvectionFD : public GenericLabOperator
{
private:
	double dt;
	
public:
	OperatorAdvectionFD(double dt) : dt(dt)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	~OperatorAdvectionFD() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double dh = info.h_gridpoint;
		const double invdh = dt*.5/dh;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				o(ix,iy).tmpU = lab(ix,iy).u + lab(ix,iy).u * invdh * (lab(ix+1,iy).u - lab(ix-1,iy).u) + lab(ix,iy).v * invdh * (lab(ix,iy+1).u - lab(ix,iy-1).u);
				o(ix,iy).tmpV = lab(ix,iy).v + lab(ix,iy).u * invdh * (lab(ix+1,iy).v - lab(ix-1,iy).v) + lab(ix,iy).v * invdh * (lab(ix,iy+1).v - lab(ix,iy-1).v);
			}
	}
};

#endif
