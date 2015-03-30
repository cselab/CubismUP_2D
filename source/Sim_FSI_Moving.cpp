//
//  Sim_FSI_Moving.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/26/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Sim_FSI_Moving.h"

#include "ProcessOperatorsOMP.h"
#include "OperatorDivergence.h"
#include "OperatorVorticity.h"
#include "PoissonSolverScalarFFTW.h"
#include "OperatorGradP.h"

#include "CoordinatorIC.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPenalization.h"
#include "CoordinatorComputeShape.h"
#include "CoordinatorPressure.h"

void Sim_FSI_Moving::_diagnostics()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	double drag = 0;
	double volume = 0;
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for schedule(static) reduction(+:drag) reduction(+:volume)
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				if (b(ix,iy).chi>0)
				{
					drag += (b(ix,iy).u-uBody[0]) * b(ix,iy).chi;
					volume += b(ix,iy).chi;
				}
			}
	}
	
	drag *= dh*dh*lambda;
	volume *= dh*dh;
	
	const double cD = 2*drag/(uBody[0]*uBody[0]*shape->getCharLength());
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	if (verbose)
		cout << step << " " << _nonDimensionalTime() << " " << bpdx << " " << dt << " " << dtCFL << " " << dtFourier << " " << drag << " " << lambda << endl;
	myfile << step << " " << _nonDimensionalTime() << " " << bpdx << " " << dt << " " << dtCFL << " " << dtFourier << " " << cD << " " << lambda << endl;
}

void Sim_FSI_Moving::_ic()
{
	CoordinatorIC coordIC(shape,0,grid);
	profiler.push_start(coordIC.getName());
	coordIC(0);
	
	stringstream ss;
	ss << path2file << "-IC.vti";
	dumper.Write(*grid, ss.str());
	profiler.pop_stop();
}

double Sim_FSI_Moving::_nonDimensionalTime()
{
	return time*abs(uBody[0])/shape->getCharLength();
}

Sim_FSI_Moving::Sim_FSI_Moving(const int argc, const char ** argv) : Simulation_FSI(argc, argv), uBody{0,0}, omegaBody(0), re(0), nu(0), dtBody(0), dtCFL(0), dtFourier(0)
{
	int rank = 0;
#ifdef _MULTIGRID_
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif // _MULTIGRID_
	
	if (rank==0)
	{
		cout << "====================================================================================================================\n";
		cout << "\t\t\tFlow past a moving cylinder\n";
		cout << "====================================================================================================================\n";
	}
	
	// simulation settings
	re = parser("-Re").asDouble(100);
	
	uBody[0] = - parser("-uBody").asDouble(0.1);
	Real center[2] = {.85,.5};
	shape->setPosition(center);
	nu = shape->getCharLength()*abs(uBody[0])/re;
	
	_ic();
	
	pipeline.clear();
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid));
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
	pipeline.push_back(new CoordinatorPressureSimple<Lab>(grid));
	pipeline.push_back(new CoordinatorPenalization(&uBody[0], &uBody[1], &omegaBody, shape, lambda, grid));
	pipeline.push_back(new CoordinatorComputeShape(&uBody[0], &uBody[1], &omegaBody, shape, grid));
	
	cout << "Coordinator/Operator ordering:\n";
	for (int c=0; c<pipeline.size(); c++)
		cout << "\t" << pipeline[c]->getName() << endl;
	
	assert(uBody[1] == 0);
}

Sim_FSI_Moving::~Sim_FSI_Moving()
{
}

void Sim_FSI_Moving::simulate()
{
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
	
	time = 0;
	double nextDumpTime = 0;
	double maxU = uBody[0];
	while (true)
	{
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		
		// choose dt (CFL, Fourier)
		profiler.push_start("DT");
		maxU = findMaxUOMP(vInfo,*grid);
		dtFourier = .1*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
		dtCFL     = .1*vInfo[0].h_gridpoint/abs(maxU);
		dtBody    = .1*vInfo[0].h_gridpoint/abs(uBody[0]);
		dt = min(min(dtCFL,dtFourier),dtBody);
		if (dumpTime>0)
			dt = min(dt,nextDumpTime-_nonDimensionalTime());
		if (endTime>0)
			dt = min(dt,endTime);
		if (verbose)
			cout << "dt (Fourier, CFL, body): " << dtFourier << " " << dtCFL << " " << dtBody << endl;
		profiler.pop_stop();
		
		for (int c=0; c<pipeline.size(); c++)
		{
			profiler.push_start(pipeline[c]->getName());
			(*pipeline[c])(dt);
			profiler.pop_stop();
		}
		
		time += dt;
		step++;
		
		
		// compute diagnostics
		if (step % 10 == 0)
		{
			profiler.push_start("Diagnostics");
			_diagnostics();
			profiler.pop_stop();
		}
		
		//dump some time steps every now and then
		profiler.push_start("Dump");
		_dump(nextDumpTime);
		profiler.pop_stop();
		
		
		if(step % 1000 == 0)
			profiler.printSummary();
		
		// check nondimensional time
		if ((endTime>0 && abs(_nonDimensionalTime()-endTime) < 10*std::numeric_limits<Real>::epsilon()) || (nsteps!=0 && step>=nsteps))
		{
			profiler.push_start("Dump");
			stringstream ss;
			ss << path2file << "-Final.vti";
			cout << ss.str() << endl;
			
			dumper.Write(*grid, ss.str());
			
			Layer vorticity(sizeX,sizeY,1);
			processOMP<Lab, OperatorVorticity>(vorticity,vInfo,*grid);
			stringstream sVort;
			sVort << path2file << "Vorticity-Final.vti";
			dumpLayer2VTK(step,sVort.str(),vorticity,1);
			profiler.pop_stop();
			
			profiler.printSummary();
			
			exit(0);
		}
	}
}