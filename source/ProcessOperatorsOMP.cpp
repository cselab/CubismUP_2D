//
//  ProcessOperatorsOMP.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/9/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "ProcessOperatorsOMP.h"
#include <cmath>

void updateOMP(vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	const int N = myInfo.size();
	
#pragma omp parallel for schedule(static)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				b(ix,iy).u = b(ix,iy).tmpU;
				b(ix,iy).v = b(ix,iy).tmpV;
#ifdef _MULTIPHASE_
				b(ix,iy).rho = b(ix,iy).tmp;
#endif // _MULTIPHASE_
			}
	}
};

void updateRhoOMP(vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	const int N = myInfo.size();
	
#pragma omp parallel for schedule(static)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				b(ix,iy).rho = b(ix,iy).tmp;
	}
};

void updatePressuresOMP(vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	const int N = myInfo.size();
	
#pragma omp parallel for schedule(static)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				b(ix,iy).pOld = b(ix,iy).p;
				b(ix,iy).p    = b(ix,iy).divU;
			}
	}
};

void resetOMP(vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	const int N = myInfo.size();
	
#pragma omp parallel for schedule(static)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				b(ix,iy).tmpU = 0;
				b(ix,iy).tmpV = 0;
#ifdef _MULTIPHASE_
				b(ix,iy).tmp = 0;
#endif // _MULTIPHASE_
			}
	}
};

void resetRhoOMP(vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	const int N = myInfo.size();
	
#pragma omp parallel for schedule(static)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				b(ix,iy).tmp = 0;
	}
};

double findMaxUOMP(vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	double maxU = 0;
	const int N = myInfo.size();
	
#pragma omp parallel for schedule(static) reduction(max:maxU)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				maxU = max(maxU,(double)abs(b(ix,iy).u));
				maxU = max(maxU,(double)abs(b(ix,iy).v));
			}
	}
	
	return maxU;
};

void computeBodyVelocity(vector<BlockInfo>& myInfo, FluidGrid & grid, Real ub[2], Real& angularU, Real rhoS, Real g[2], Real dt, Real lambda)
{
	double centerTmpX = 0;
	double centerTmpY = 0;
	double mass = 0;
	double u = 0;
	double v = 0;
	double momOfInertia = 0;
	double angularMomentum = 0;
	const int N = myInfo.size();
	
#pragma omp parallel for schedule(static) reduction(+:centerTmpX) reduction(+:centerTmpY) reduction(+:mass)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
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
		BlockInfo info = myInfo[i];
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
	
	ub[0] = u / mass;
	ub[1] = v / mass;
	angularU = angularMomentum / momOfInertia;
	
	/*/
#pragma omp parallel for schedule(static) reduction(+:u) reduction(+:v)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
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
};

void computeForcesFromVorticity(vector<BlockInfo>& myInfo, FluidGrid & grid, Real ub[2], Real oldAccVort[2], Real rhoS)
{
	Real mU = 0;
	Real mV = 0;
	Real mass = 0;
	const int N = myInfo.size();
	
#pragma omp parallel for schedule(static) reduction(+:mU) reduction(+:mV) reduction(+:mass)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		Real h = info.h_gridpoint;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				Real p[2];
				info.pos(p, ix, iy);
				
				mU += - (1-b(ix,iy).chi) * p[1] * b(ix,iy).tmp * b(ix,iy).rho;
				mV -= - (1-b(ix,iy).chi) * p[0] * b(ix,iy).tmp * b(ix,iy).rho;
				mass += b(ix,iy).chi * rhoS;
			}
	}
	
	ub[0] += (mU-oldAccVort[0]) / mass;
	ub[1] += (mV-oldAccVort[1]) / mass;
	oldAccVort[0] = mU;
	oldAccVort[1] = mV;
}
