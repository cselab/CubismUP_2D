//
//  TestPressure.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/9/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestPressure.h"
#include "OperatorDivergence.h"
#include "OperatorGradP.h"
#include "PoissonSolverScalarFFTW.h"
#include "ProcessOperatorsOMP.h"
#include "LayerToVTK.h"
#include <sstream>
#include <cmath>
#ifdef _MULTIGRID_
#include "MultigridHypre.h"
#endif

void initializeNUMA(Layer& layer)
{
	for(int ic = 0; ic < layer.nDim; ic++)
#pragma omp parallel for
		for(int iy = 0; iy < layer.sizeY; iy++)
		{
			for(int ix=0; ix<layer.sizeX; ix++)
			{
				layer.data[ic*layer.sizeX*layer.sizeY + iy*layer.sizeX + ix] = 0;
				
			}
		}
}

class BS4
{
public:
	static inline Real eval(Real x)
	{
		const Real t = fabs(x);
		
		if (t>2) return 0;
		
		if (t>1) return pow(2-t,3)/6;
		
		return (1 + 3*(1-t)*(1 + (1-t)*(1 - (1-t))))/6;
	}
};

void TestPressure::_ic()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					double p[3];
					info.pos(p, ix, iy);
					p[0] = p[0]*2.-1.;
					p[1] = p[1]*2.-1.;
					
					if (ic==0)
					{
						double x = p[0]*M_PI;
						double y = p[1]*M_PI;
						
						b(ix, iy).u   = 1./(4.*M_PI*M_PI)*cos(x); // expected solution
						b(ix, iy).divU = -cos(x); // rhs
					}
					else if (ic==1)
					{
						const Real IrI  = sqrt(p[0]*p[0] + p[1]*p[1]);
						const double strength = 100./(1+IrI*IrI)*BS4::eval(IrI/0.5*2.5);
						
						b(ix, iy).rho = 1./(4.*M_PI*M_PI)*cos(p[0]*M_PI);
						b(ix, iy).u   = -p[1]*strength;
						b(ix, iy).v   =  p[0]*strength;
						b(ix, iy).chi = 0;
						
						b(ix, iy).divU = -cos(p[0]*M_PI);
					}
					else if (ic==2)
					{
						// this is a test for:
						//	0-dirichlet on x=0
						//	0-neumann on x=1
						
						const int size = 1+1/dh;
						const int by = info.index[1]*FluidBlock::sizeY;
						p[1] = (by+iy+1.5)/(double)size;
						double y = p[1]*M_PI_2;
						b(ix,iy).rho = 1;
						b(ix,iy).divU = -M_PI*M_PI_4 * cos(y);
						b(ix,iy).u = cos(y);
					}
				}
	}
	
	stringstream ss;
	ss << path2file << "-ic" << ic << "-bpd" << bpd << "-IC.vti";
	
	dumper.Write(*grid, ss.str());
}

TestPressure::TestPressure(const int argc, const char ** argv, const int solver, const int ic, const int bpd) : Test(argc, argv), solver(solver), ic(ic), bpd(bpd)
{
#ifdef _MULTIGRID_
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
#endif
	
	grid = new FluidGrid(bpd,bpd,1);
	
	// output settings
	path2file = parser("-file").asString("../data/testPressure");
	
	// setup initial condition
#ifdef _MULTIGRID_
	if (rank==0)
#endif
		_ic();
}

TestPressure::~TestPressure()
{
	delete grid;
}

void TestPressure::run()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const int size = bpd * FluidBlock::sizeX;
	
	if (ic==1 && rank==0)
		processOMP<Lab, OperatorDivergence>(1, vInfo, *grid);
	
	if (solver==0)
	{
#ifdef _SP_COMP_
		PoissonSolverScalarFFTW<FluidGrid, StreamerDiv> pressureSolver(NTHREADS);
		pressureSolver.solve(*grid,false);
#else
		cout << "FFTW double precision not supported - aborting now!\n";
		abort();
#endif
	}
	else if (solver==1)
	{
#ifdef _SP_COMP_
		PoissonSolverScalarFFTW<FluidGrid, StreamerDiv> pressureSolver(NTHREADS);
		pressureSolver.solve(*grid,true);
#else
		cout << "FFTW double precision not supported - aborting now!\n";
		abort();
#endif
	}
#ifdef _MULTIGRID_
	else if (solver==2)
	{
		const bool bConstantCoefficients = true;
		MultigridHypre mg;
		mg.setup(grid, bConstantCoefficients, rank, nprocs);
		mg();
	}
	else if (solver==3)
	{
		cout << "Don't use it for this tests\n";
		abort();
		const bool bConstantCoefficients = false;
		MultigridHypre mg;
		mg.setup(grid, bConstantCoefficients, rank, nprocs);
		mg();
	}
#endif
	
	if (ic==1 && rank==0)
		processOMP<Lab, OperatorGradP>(1, vInfo, *grid);
	
	if (rank==0)
	{
		stringstream ss;
		ss << path2file << "-solver" << solver << "-ic" << ic << "-bpd" << bpd << ".vti";
	
		dumper.Write(*grid, ss.str());
	}
}

void TestPressure::check()
{
#ifdef _MULTIGRID_
	if (rank==0)
#endif
	{
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		const int size = bpd * FluidBlock::sizeX;
		
		Layer divergence(size,size,1);
		processOMP<Lab, OperatorDivergenceLayer>(divergence,vInfo,*grid);
		
		
		cout << "\tErrors (Linf, L1, L2):\t";
		double Linf = 0.;
		double L1 = 0.;
		double L2 = 0.;
		
		stringstream ss;
		ss << path2file << "_diagnostics.dat";
		ofstream myfile(ss.str(), fstream::app);
		
		const double dh = vInfo[0].h_gridpoint;
		
#pragma omp parallel for reduction(max:Linf) reduction(+:L1) reduction(+:L2)
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					double error;
					if (ic==0)
						error = b(ix,iy).divU - b(ix,iy).u;
					else
						error = divergence(ix,iy);
					
					Linf = max(Linf,abs(error));
					L1 += abs(error);
					L2 += error*error;
				}
		}
		
		L1 *= dh*dh;
		L2 = sqrt(L2)*dh;
		cout << "\t" << Linf << "\t" << L1 << "\t" << L2 << endl;
		myfile << size << " " << Linf << " " << L1 << " " << L2 << endl;
	}
}
