//
//  OperatorVarCoeffPoisson.h
//  CubismUP_2D
//
//	Operates on
//		tmpU, tmp, divU
//
//	Using Jacobi iterative method
//
//  Created by Christian Conti on 1/23/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_OperatorVarCoeffPoisson_h
#define CubismUP_2D_OperatorVarCoeffPoisson_h

#include <cmath>
#include "GenericOperator.h"

// this does a simple iteration
class OperatorVarCoeffPoisson : public GenericLabOperator
{
private:
	double relaxParam;
	
	inline void _mean(Real c, Real e, Real w, Real n, Real s, Real& avgE, Real& avgW, Real& avgN, Real& avgS)
	{
		avgE = .5 * (c + e);
		avgW = .5 * (c + w);
		avgN = .5 * (c + n);
		avgS = .5 * (c + s);
	}
	
	inline void _harmonicAvg(Real c, Real e, Real w, Real n, Real s, Real& avgE, Real& avgW, Real& avgN, Real& avgS)
	{
		avgE = 2. * c * e / (c + e);
		avgW = 2. * c * w / (c + w);
		avgN = 2. * c * n / (c + n);
		avgS = 2. * c * s / (c + s);
	}
	
public:
	OperatorVarCoeffPoisson(double relaxParam=1) : relaxParam(relaxParam)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	
	~OperatorVarCoeffPoisson() {}
	
	// this implementation is still biased: some blocklabs might not be loaded yet when another block is finished (propagation of information)
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo & info, BlockType & o)
	{
		const Real dh2 = info.h_gridpoint * info.h_gridpoint;
		
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				FluidElement& phi = lab(ix,iy);
				FluidElement& phiN = lab(ix,iy+1);
				FluidElement& phiS = lab(ix,iy-1);
				FluidElement& phiE = lab(ix+1,iy);
				FluidElement& phiW = lab(ix-1,iy);
				
				Real rhoE, rhoW, rhoN, rhoS;
				//_mean(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
				_harmonicAvg(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
				
				Real factor = 1./rhoE + 1./rhoW + 1./rhoN + 1./rhoS;
				o(ix,iy).tmpU = 1./factor * (phiE.tmp/rhoE + phiW.tmp/rhoW + phiN.tmp/rhoN + phiS.tmp/rhoS - dh2*phi.divU);
				//Real factor = rhoE + rhoW + rhoN + rhoS;
				//o(ix,iy).tmpU = (1-relaxParam)*phi.tmp + relaxParam/factor * (phiE.tmp*rhoE + phiW.tmp*rhoW + phiN.tmp*rhoN + phiS.tmp*rhoS - dh2*phi.divU);
				//o(ix,iy).tmpU = .25 * (phiE.tmp + phiW.tmp + phiN.tmp + phiS.tmp - dh2*phi.divU);
			}
	}
};

template<typename Lab, typename Kernel>
void processOMP_Jacobi(vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	BlockInfo * ary = &myInfo.front();
	const int N = myInfo.size();
	const int niter = 100000;
	int iter = 0;
	const double tol = 1e-7;
	double Linf = 0.;
	double L1 = 0.;
	double L2 = 0.;
	
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
	// 0. initialize guess - take the existing one (IC for time 0, previous solution for time>0)
	
	while (iter<niter)
	{
		// 1. iterate
#pragma omp parallel
		{
			Kernel kernel;
			
			Lab mylab;
			mylab.prepare(grid, kernel.stencil_start, kernel.stencil_end, true);
			
#pragma omp for schedule(static)
			for(int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0, false);
				
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
		/*
		if (iter%100==0)
		{
		stringstream ssPre;
		ssPre << "../data/PoissonPre" << iter << ".vti" ;
		dumper.Write(grid, ssPre.str());
		}
		 */
		
		// 2. check convergence
		Linf = 0;
		L1 = 0;
		L2 = 0;
#pragma omp parallel for schedule(static) reduction(max:Linf) reduction(+:L1) reduction(+:L2)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = myInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					double diff = b(ix,iy).tmp - b(ix,iy).tmpU;
					
					Linf = max(Linf,abs(diff));
					L1 += abs(diff);
					L2 += diff*diff;
				}
		}
		const double dh = myInfo[0].h_gridpoint;
		L1 *= dh*dh;
		L2 = sqrt(L2)*dh;
		//cout << Linf << " " << L1 << " " << L2 << endl;
		
		// 3. overwrite old results with new ones
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = myInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					b(ix,iy).tmp = b(ix,iy).tmpU;
				}
		}
		
		if (Linf < tol)
		break;
		
		iter++;
	}
	cout << "Converged in " << iter << " iterations (Linf is " << Linf << ")\n";
	
	// put final results into divU
#pragma omp parallel for schedule(static)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix)
		{
			b(ix,iy).divU = b(ix,iy).tmp;
		}
	}
}

template<typename Lab, typename Kernel>
void processOMP_SRJ(vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	BlockInfo * ary = &myInfo.front();
	const int N = myInfo.size();
	const int niter = 100000;
	int iter = 0;
	const double tol = 1e-8;
	double Linf = 0.;
	double L1 = 0.;
	double L2 = 0.;
	
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
	/*
	const int P = 1;
	const int M = 1;
	const double relaxParam[1] = { 1 };
	const int Q[1] = { 1 };
	/*/
	const int P = 2;
	//const int M = 2;
	//const double relaxParam[2] = { 3.414213, 0.585786 };
	//const int Q[2] = { 1, 1 };
	//const int M = 64;
	//const double relaxParam[2] = { 190.2, 0.9532 };
	//const int Q[2] = { 1, 63 };
	const int M = 16;
	const double relaxParam[2] = { 32.60, 0.8630 };
	const int Q[2] = { 1, 15 };
	/*
	 // these are relaxation parameters for grids of at least 64^2
	const int P = 5;
	const int M = 158;
	const double relaxParam[5] = { 1228.8, 220.14, 26.168, 3.1668, 0.63890 };
	const int Q[5] = { 1, 3, 10, 38, 106 };
	 */
	
	// 0. initialize guess - take the existing one (IC for time 0, previous solution for time>0)
	
	while (iter*M<niter)
	{
		for (int level=0; level<P; level++)
		{
			for (int mIter = 0; mIter < Q[level]; mIter++)
			{
				// 1. iterate
#pragma omp parallel
				{
					Kernel kernel(relaxParam[level]);
					
					Lab mylab;
					mylab.prepare(grid, kernel.stencil_start, kernel.stencil_end, true);
					
#pragma omp for schedule(static)
					for(int i=0; i<N; i++)
					{
						mylab.load(ary[i], 0, false);
						
						kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
					}
				}
			}
		}
		
		/*
		{
			stringstream ssPre;
			ssPre << "../data/PoissonPre" << iter << ".vti" ;
			dumper.Write(grid, ssPre.str());
		}
		 */
		
		// 2. check convergence
		Linf = 0;
		L1 = 0;
		L2 = 0;
#pragma omp parallel for schedule(static) reduction(max:Linf) reduction(+:L1) reduction(+:L2)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = myInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					double diff = b(ix,iy).tmp - b(ix,iy).tmpU;
					
					Linf = max(Linf,abs(diff));
					L1 += abs(diff);
					L2 += diff*diff;
				}
		}
		const double dh = myInfo[0].h_gridpoint;
		L1 *= dh*dh;
		L2 = sqrt(L2)*dh;
		
		if (Linf < tol)
			break;
		
		// 3. overwrite old results with new ones
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = myInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					b(ix,iy).tmp = b(ix,iy).tmpU;
				}
		}
		
		iter++;
	}
	cout << "Converged in " << iter*M << " iterations (Linf is " << Linf << ")\n";
}

#endif
