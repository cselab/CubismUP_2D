//
//  Sim_FSI_Fixed.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Sim_FSI_Fixed.h"

#include "ProcessOperatorsOMP.h"
#include "OperatorDivergence.h"
#include "OperatorVorticity.h"
#include "OperatorGradP.h"
#include "OperatorComputeShape.h"

#include "CoordinatorIC.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPenalization.h"
#include "CoordinatorPressure.h"


void Sim_FSI_Fixed::_diagnostics()
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
					drag += b(ix,iy).u * b(ix,iy).chi;
					volume += b(ix,iy).chi;
				}
			}
	}
	
	drag *= dh*dh*lambda;
	volume *= dh*dh;
	
	const double cD = 2.*drag/(uinf*uinf*shape->getCharLength());
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	if (verbose)
		cout << step << " " << _nonDimensionalTime() << " " << bpdx << " " << dt << " " << dtCFL << " " << dtFourier << " " << drag << " " << lambda << endl;
	myfile << step << " " << _nonDimensionalTime() << " " << bpdx << " " << dt << " " << dtCFL << " " << dtFourier << " " << cD << " " << lambda << endl;
}

void Sim_FSI_Fixed::_ic()
{
	CoordinatorIC coordIC(shape,uinf,grid);
	profiler.push_start(coordIC.getName());
	coordIC(0);
	
	stringstream ss;
	ss << path2file << "-IC.vti";
	dumper.Write(*grid, ss.str());
	profiler.pop_stop();
}

double Sim_FSI_Fixed::_nonDimensionalTime()
{
	return time*abs(uinf)/shape->getCharLength();
}

Sim_FSI_Fixed::Sim_FSI_Fixed(const int argc, const char ** argv) : Simulation_FSI(argc, argv), uinf(0), re(0), nu(0), dtCFL(0), dtFourier(0)
{
	int rank = 0;
#ifdef _MULTIGRID_
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif // _MULTIGRID_
	
	if (rank==0)
	{
		cout << "====================================================================================================================\n";
		cout << "\t\t\tFlow past a fixed body\n";
		cout << "====================================================================================================================\n";
	}
	
	// simulation settings
	uinf = parser("-uinf").asDouble(0.1);
	re = parser("-Re").asDouble(100);
	
	Real center[2] = {.15,.5};
	shape->setPosition(center);
	nu = shape->getCharLength()*uinf/re;
	
	_ic();
	
	pipeline.clear();
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid));
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
	pipeline.push_back(new CoordinatorPressureSimple<Lab>(grid));
	pipeline.push_back(new CoordinatorPenalizationFixed(shape, lambda, grid));
	
	cout << "Coordinator/Operator ordering:\n";
	for (int c=0; c<pipeline.size(); c++)
		cout << "\t" << pipeline[c]->getName() << endl;
}

Sim_FSI_Fixed::~Sim_FSI_Fixed()
{
}

void Sim_FSI_Fixed::simulate()
{
	const Real uBody[2] = {0,0};
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
	
	time = 0;
	double nextDumpTime = 0;
	double maxU = uinf;
    while (true)
	{
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		
		// choose dt (CFL, Fourier)
		profiler.push_start("DT");
		maxU = findMaxUOMP(vInfo,*grid);
		dtFourier = .1*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
		dtCFL     = .1*vInfo[0].h_gridpoint/abs(maxU);
		dt = min(dtCFL,dtFourier);
		if (dumpTime>0)
			dt = min(dt,nextDumpTime-_nonDimensionalTime());
		if (endTime>0)
			dt = min(dt,endTime-_nonDimensionalTime());
		if (verbose)
			cout << "dt (Fourier, CFL): " << dtFourier << " " << dtCFL << endl;
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
		
			
		if (step % 1000 == 0)
			profiler.printSummary();
		
		// check nondimensional time
		if ((endTime>0 && abs(_nonDimensionalTime()-endTime) < 10*std::numeric_limits<Real>::epsilon()) || (nsteps!=0 && step>=nsteps))
		{
			if (verbose)
				cout << "Finished at time " << _nonDimensionalTime() << " (target end time " << endTime << ") in " << step << " step of " << nsteps << endl;
			
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