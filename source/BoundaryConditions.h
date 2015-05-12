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
	
	template<int dir, int side>
	void _setup()
	{
		s[0] =	dir==0? (side==0? stencilStart[0]: TBlock::sizeX) : stencilStart[0];
		s[1] =	dir==1? (side==0? stencilStart[1]: TBlock::sizeY) : stencilStart[1];
		s[2] =	dir==2? (side==0? stencilStart[2]: TBlock::sizeZ) : stencilStart[2];
		
		e[0] =	dir==0? (side==0? 0: TBlock::sizeX + stencilEnd[0]-1) : TBlock::sizeX +  stencilEnd[0]-1;
		e[1] =	dir==1? (side==0? 0: TBlock::sizeY + stencilEnd[1]-1) : TBlock::sizeY +  stencilEnd[1]-1;
		e[2] =	dir==2? (side==0? 0: TBlock::sizeZ + stencilEnd[2]-1) : TBlock::sizeZ +  stencilEnd[2]-1;
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
				(*this)(ix,iy).rk2u = p.rk2u;
				(*this)(ix,iy).rk2v = p.rk2v;
			}
	}
	
	template<int dir, int side>
	void applyBC_neumann()
	{
		_setup<dir,side>();
		
		for(int iy=s[1]; iy<e[1]; iy++)
			for(int ix=s[0]; ix<e[0]; ix++)
			{
				(*this)(ix,iy) = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1):ix,
										 dir==1? (side==0? 0:TBlock::sizeY-1):iy);
			}
	}
	
	template<int dir, int side>
	void applyBC_mixedBottom()
	{
		_setup<dir,side>();
		
		for(int iy=s[1]; iy<e[1]; iy++)
			for(int ix=s[0]; ix<e[0]; ix++)
			{
				(*this)(ix,iy).rho = 1;
				(*this)(ix,iy).chi = 0;
				(*this)(ix,iy).u   = 0;
				(*this)(ix,iy).v   = 0;
				(*this)(ix,iy).rk2u = 0;
				(*this)(ix,iy).rk2v = 0;
				(*this)(ix,iy).p   = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1-ix+s[0]):ix,
											 dir==1? (side==0? 0:TBlock::sizeY-1-iy+s[1]):iy).p;
				(*this)(ix,iy).pOld = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1-ix+s[0]):ix,
											  dir==1? (side==0? 0:TBlock::sizeY-1-iy+s[1]):iy).pOld;
				(*this)(ix,iy).divU = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1-ix+s[0]):ix,
											  dir==1? (side==0? 0:TBlock::sizeY-1-iy+s[1]):iy).divU;
			}
	}
	
	template<int dir, int side>
	void applyBC_mixedTop()
	{
		_setup<dir,side>();
		
		for(int iy=s[1]; iy<e[1]; iy++)
			for(int ix=s[0]; ix<e[0]; ix++)
			{
				(*this)(ix,iy).rho = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1-ix+s[0]):ix,
											 dir==1? (side==0? 0:TBlock::sizeY-1-iy+s[1]):iy).rho;
				(*this)(ix,iy).chi = 0;
				(*this)(ix,iy).u   = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1-ix+s[0]):ix,
											 dir==1? (side==0? 0:TBlock::sizeY-1-iy+s[1]):iy).u;
				(*this)(ix,iy).v   = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1-ix+s[0]):ix,
											 dir==1? (side==0? 0:TBlock::sizeY-1-iy+s[1]):iy).v;
				(*this)(ix,iy).rk2u = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1-ix+s[0]):ix,
											 dir==1? (side==0? 0:TBlock::sizeY-1-iy+s[1]):iy).rk2u;
				(*this)(ix,iy).rk2v = (*this)(dir==0? (side==0? 0:TBlock::sizeX-1-ix+s[0]):ix,
											 dir==1? (side==0? 0:TBlock::sizeY-1-iy+s[1]):iy).rk2v;
				(*this)(ix,iy).p   = 0;
				(*this)(ix,iy).pOld = 0;
				(*this)(ix,iy).divU = 0;
			}
	}
};
