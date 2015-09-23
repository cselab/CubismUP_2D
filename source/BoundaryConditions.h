/*
 *  BoundaryConditions.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/16/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "BlockLab.h"
#include <Matrix3D.h>

template<typename TBlock, typename TElement, template<typename X> class allocator=std::allocator>
class BoundaryCondition
{
protected:
	
	int s[3], e[3];
	int stencilStart[3], stencilEnd[3];
	Matrix3D<TElement, true, allocator> * cacheBlock;
	
	// this guy includes corners!
	// what happens if they are not allocated?
	template<int dir, int side>
	void _setup()
	{
		s[0] =	dir==0 ? (side==0 ? stencilStart[0]: TBlock::sizeX) : stencilStart[0];
		s[1] =	dir==1 ? (side==0 ? stencilStart[1]: TBlock::sizeY) : stencilStart[1];
		s[2] =	dir==2 ? (side==0 ? stencilStart[2]: TBlock::sizeZ) : stencilStart[2];
		
		e[0] =	dir==0 ? (side==0 ? 0: TBlock::sizeX + stencilEnd[0]-1) : TBlock::sizeX +  stencilEnd[0]-1;
		e[1] =	dir==1 ? (side==0 ? 0: TBlock::sizeY + stencilEnd[1]-1) : TBlock::sizeY +  stencilEnd[1]-1;
		e[2] =	dir==2 ? (side==0 ? 0: TBlock::sizeZ + stencilEnd[2]-1) : TBlock::sizeZ +  stencilEnd[2]-1;
	}
public:
	
	BoundaryCondition(const int ss[3], const int se[3], Matrix3D<TElement, true, allocator> * cacheBlock):
	cacheBlock(cacheBlock)
	{
		s[0]=s[1]=s[2]=0;
		e[0]=e[1]=e[2]=0;
		
		stencilStart[0] = ss[0];
		stencilStart[1] = ss[1];
		stencilStart[2] = ss[2];
		
		stencilEnd[0] = se[0];
		stencilEnd[1] = se[1];
		stencilEnd[2] = se[2];
	}
	
	TElement& operator()(int ix, int iy)
	{
		return cacheBlock->Access(ix-stencilStart[0],iy-stencilStart[1], 0);
	}
	
	template<int dir, int side>
	void applyBC_dirichlet(const TElement& p)
	{
		_setup<dir,side>();
		
		for(int iy=s[1]; iy<e[1]; iy++)
			for(int ix=s[0]; ix<e[0]; ix++)
			{
				(*this)(ix,iy).rho = p.rho;
				(*this)(ix,iy).u   = p.u;
				(*this)(ix,iy).v   = p.v;
				(*this)(ix,iy).chi = p.chi;
				(*this)(ix,iy).tmp = p.tmp;
			}
	}
	
	template<int dir, int side>
	void applyBC_neumann()
	{
		_setup<dir,side>();
        abort();
		for(int iy=s[1]; iy<e[1]; iy++)
			for(int ix=s[0]; ix<e[0]; ix++)
            {
                // the Neumann BC is incorrect!
				(*this)(ix,iy) = (*this)(dir==0 ? (side==0 ? 0 : TBlock::sizeX-1) : ix,
										 dir==1 ? (side==0 ? 0 : TBlock::sizeY-1) : iy);
			}
	}
	
	template<int dir, int side>
	void applyBC_mixedBottom(const TElement& p)
	{
		assert(dir==1);
        assert(side==0);
		
		_setup<dir,side>();
		
		for(int iy=s[1]; iy<e[1]; iy++)
			for(int ix=s[0]; ix<e[0]; ix++)
            {
                // the Neumann BC is incorrect! This might explain the issues with the high order diffusion!
                (*this)(ix,iy).rho  = (*this)(ix, -iy).rho;
				(*this)(ix,iy).chi  = p.chi;
				(*this)(ix,iy).u    = p.u;
				(*this)(ix,iy).v    = -p.v;
				(*this)(ix,iy).p    = (*this)(ix, -iy).p;
				(*this)(ix,iy).pOld = (*this)(ix, -iy).pOld;
				(*this)(ix,iy).divU = (*this)(ix, -iy).divU;
			}
	}
	
	template<int dir, int side>
	void applyBC_mixedTop(const TElement& p)
	{
		assert(dir==1);
        assert(side==1);
		
		_setup<dir,side>();
		
		for(int iy=s[1]; iy<e[1]; iy++)
			for(int ix=s[0]; ix<e[0]; ix++)
			{
                // are the Neumann BC correct?
				(*this)(ix,iy).rho  = (*this)(ix, TBlock::sizeY-2-iy+s[1]).rho;
				(*this)(ix,iy).chi  = p.chi;
				(*this)(ix,iy).u    = (*this)(ix, TBlock::sizeY-2-iy+s[1]).u;
				(*this)(ix,iy).v    = (*this)(ix, TBlock::sizeY-2-iy+s[1]).v;
				(*this)(ix,iy).p    = p.p;
				(*this)(ix,iy).pOld = p.pOld;
				(*this)(ix,iy).divU = p.divU;
			}
	}
	
	template<int dir, int side>
	void applyBC_vortex(const BlockInfo info)
	{
		_setup<dir,side>();
		
		for(int iy=s[1]; iy<e[1]; iy++)
			for(int ix=s[0]; ix<e[0]; ix++)
			{
				/*
				(*this)(ix,iy).rho = 1;
				(*this)(ix,iy).u = 0;
				(*this)(ix,iy).v = 0;
				(*this)(ix,iy).p = 0;
				*/
				 //*
				double p[3];
				info.pos(p, ix, iy);
				
				p[0] = p[0]*2.-1.;
				p[1] = p[1]*2.-1.;
                
                const Real r = sqrt(p[0]*p[0] + p[1]*p[1]);
				const Real invR = 1./r;
                
                (*this)(ix,iy).rho = r;
                (*this)(ix,iy).u   =   sin(p[1])*cos(r*M_PI/2)*invR;//-p[1];//
                (*this)(ix,iy).v   =  -sin(p[0])*cos(r*M_PI/2)*invR;// p[0];//
				(*this)(ix,iy).chi = 0;
				 // what about pressure?
				 //*/
			}
	}
};
