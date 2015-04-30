/*
 *  BlockInfo.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 5/24/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <cstdlib>

struct BlockInfo
{
	long long blockID;
	void * ptrBlock;
	int index[3];

	double origin[3];
	double h, h_gridpoint;
	
	template <typename T>
	inline void pos(T p[2], int ix, int iy) const
	{
#ifndef _VERTEXCENTERED_
		p[0] = origin[0] + h_gridpoint*(ix+0.5);
		p[1] = origin[1] + h_gridpoint*(iy+0.5);
#else
		p[0] = origin[0] + h_gridpoint*ix;
		p[1] = origin[1] + h_gridpoint*iy;
#endif
	}
	
	template <typename T>
	inline void pos(T p[3], int ix, int iy, int iz) const
	{
#ifndef _VERTEXCENTERED_
		p[0] = origin[0] + h_gridpoint*(ix+0.5);
		p[1] = origin[1] + h_gridpoint*(iy+0.5);
		p[2] = origin[2] + h_gridpoint*(iz+0.5);
#else
		p[0] = origin[0] + h_gridpoint*ix;
		p[1] = origin[1] + h_gridpoint*iy;
		p[2] = origin[2] + h_gridpoint*iz;
#endif
	}
	
	BlockInfo(long long ID, const int idx[3], const double pos[3], const double spacing, double h_gridpoint_, void * ptr=NULL):
	blockID(ID), ptrBlock(ptr)
	{
		index[0] = idx[0];
		index[1] = idx[1];
		index[2] = idx[2];
		
		origin[0] = pos[0];
		origin[1] = pos[1];
		origin[2] = pos[2];
		
		h = spacing;
		h_gridpoint = h_gridpoint_;
	}
	
	BlockInfo():blockID(-1), ptrBlock(NULL) {}
};
