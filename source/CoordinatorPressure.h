//
//  CoordinatorPressure.h
//  CubismUP_2D
//
//  Created by Christian Conti on 3/30/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_CoordinatorPressure_h
#define CubismUP_2D_CoordinatorPressure_h

#include "GenericCoordinator.h"
#include "OperatorDivergence.h"
#include "OperatorGradP.h"
#include "PoissonSolverScalarFFTW.h"
#ifdef _MULTIGRID_
#include "MultigridHypre.h"
#endif // _MULTIGRID_

template <typename Lab>
class CoordinatorPressure : public GenericCoordinator
{
protected:
	const int rank, nprocs;
	const double minRho;
	int * step;
	const bool bSplit;
	
#ifdef _SPLIT_
#ifdef _SP_COMP_
	PoissonSolverScalarFFTW<FluidGrid, StreamerDiv> pressureSolver;
#endif // _SP_COMP_
#endif // _SPLIT_
	
#ifdef _MULTIGRID_
	MultigridHypre mg;
#endif // _MULTIGRID_
	
	inline void updatePressure()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					b(ix,iy).pOld = b(ix,iy).p;
					b(ix,iy).p    = b(ix,iy).divU;
				}
		}
	}
	
	template <typename Operator>
	void computeSplit(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Operator kernel(dt, minRho, *step);
			
			Lab mylab;
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
	template <typename Operator>
	void compute(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Operator kernel(dt);
			
			Lab mylab;
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
public:
	CoordinatorPressure(const double minRho, int * step, const bool bSplit, FluidGrid * grid, const int rank, const int nprocs) : GenericCoordinator(grid), rank(rank), nprocs(nprocs), minRho(minRho), step(step), bSplit(bSplit)
#ifdef _SPLIT_
#ifdef _SP_COMP_
	, pressureSolver(NTHREADS)
#endif // _SP_COMP_
#endif // _SPLIT_
	{
	}
	
	void operator()(const double dt)
	{
		// need an interface that is the same for all solvers - this way the defines can be removed more cleanly
#ifdef _MULTIGRID_
		MPI_Barrier(MPI_COMM_WORLD);
#endif // _MULTIGRID_
		
		// pressure
#ifdef _SPLIT_
#ifdef _SP_COMP_
		computeSplit<OperatorDivergenceSplit>(dt);
		pressureSolver.solve(*grid,false);
		computeSplit<OperatorGradPSplit>(dt);
#else // _SP_COMP_
		cout << "FFTW double precision not supported - aborting now!\n";
		abort();
#endif // _SP_COMP_
#endif // _SPLIT_
#ifdef _MULTIGRID_
		if (rank==0)
			if (bSplit)
				computeSplit<OperatorDivergenceSplit>(dt);
			else
				compute<OperatorDivergence>(dt);
		mg.setup(grid, bSplit, rank, nprocs);
		mg();
		if (rank==0)
			if (bSplit)
				computeSplit<OperatorGradPSplit>(dt);
			else
				compute<OperatorGradP>(dt);
#endif // _MULTIGRID_
		
#ifdef _MULTIGRID_
		if (rank==0)
#endif // _MULTIGRID_
			updatePressure();
	}
	
	string getName()
	{
		return "Pressure";
	}
};

template <typename Lab>
class CoordinatorPressureSimple : public GenericCoordinator
{
protected:
#ifdef _SP_COMP_
	PoissonSolverScalarFFTW<FluidGrid, StreamerDiv> pressureSolver;
#endif // _SP_COMP_
	
	template <typename Operator>
	void compute(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Operator kernel(dt);
			
			Lab mylab;
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
public:
	CoordinatorPressureSimple(FluidGrid * grid) : GenericCoordinator(grid)
#ifdef _SP_COMP_
	, pressureSolver(NTHREADS)
#endif // _SP_COMP_
	{
	}
	
	void operator()(const double dt)
	{
		compute<OperatorDivergence>(dt);
#ifdef _SP_COMP_
		pressureSolver.solve(*grid,true);
#else // _SP_COMP_
		cout << "FFTW double precision not supported - aborting now!\n";
		abort();
#endif // _SP_COMP_
		compute<OperatorGradP>(dt);
	}
	
	string getName()
	{
		return "Pressure";
	}
};
#endif
