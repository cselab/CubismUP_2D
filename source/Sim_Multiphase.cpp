//
//  Sim_Multiphase.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 4/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Sim_Multiphase.h"

#include "ProcessOperatorsOMP.h"
#include "OperatorVorticity.h"

#include "CoordinatorIC.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPenalization.h"
#include "CoordinatorPressure.h"
#include "CoordinatorGravity.h"
#include "CoordinatorUpdate.h"
#include "CoordinatorCleanTmp.h"

void Sim_Multiphase::_diagnostics()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	double drag = 0;
	double volS = 0;
	double volF = 0;
	double pMin = 10;
	double pMax = 0;
	const double dh = vInfo[0].h_gridpoint;
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	if (verbose)
		cout << step << " " << time << " " << dt << " " << bpdx << endl;
	myfile << step << " " << time << " " << dt << " " << bpdx << endl;
}

void Sim_Multiphase::_dumpSettings(ostream& outStream)
{
#ifdef _MULTIGRID_
	if (rank==0)
#endif // _MULTIGRID_
	{
		outStream << "--------------------------------------------------------------------\n";
		outStream << "Physical Settings\n";
		outStream << "\tnu\t" << nu << endl;
		outStream << "\tRhoS\t" << rhoS << endl;

		outStream << "\nSimulation Settings\n";
		outStream << "\tnsteps\t" << nsteps << endl;
		outStream << "\tTend\t" << endTime << endl;
#ifdef _MULTIGRID_
		outStream << "\tsplit\t" << (bSplit ? "true" : "false") << endl;
#endif // _MULTIGRID_
		outStream << "\tpath2file\t" << path2file << endl;
		outStream << "\tbpdx\t" << bpdx << endl;
		outStream << "\tbpdy\t" << bpdy << endl;
#ifdef _PERIODIC_
		outStream << "\tBC\t\tperiodic\n";
#else // _PERIODIC_
		outStream << "\tBC\tmixed\n";
#endif // _PERIODIC_
#ifdef _MULTIGRID_
		outStream << "\tPoisson\tMultigrid\n";
#endif // _MULTIGRID_
#ifdef _SPLIT_
		outStream << "\tPoisson\tFFTW Split\n";
#endif // _SPLIT_
		outStream << "--------------------------------------------------------------------\n";
	}
}

void Sim_Multiphase::_ic()
{
#ifdef _MULTIGRID_
	if (rank==0)
#endif // _MULTIGRID_
	{
		// setup initial conditions
		/*
		Shape * shape;
		Real radius = parser("-radius").asDouble(0.1);
		Real centerOfMass[2];
		centerOfMass[0] = .5;
		if (rhoS >= 1)
			centerOfMass[1] = .85;
		else
			centerOfMass[1] = .15;
		
		bool bPeriodic[2] = {false,false};
		shape = new Disk(centerOfMass, radius, rhoS, 2, 2, bPeriodic);
		CoordinatorIC coordIC(shape,0,grid);
		profiler.push_start(coordIC.getName());
		/*/
		CoordinatorIC_RT coordIC(grid);
		profiler.push_start(coordIC.getName());
		*/
		coordIC(0);
		
		stringstream ss;
		ss << path2file << "-IC.vti";
		dumper.Write(*grid, ss.str());
		profiler.pop_stop();
		
		delete shape;
	}
}

double Sim_Multiphase::_nonDimensionalTime()
{
	return time; // how to nondimensionalize here? based on Galileo number?
}

void Sim_Multiphase::_outputSettings(ostream &outStream)
{
	outStream << "Multiphase\n";
	outStream << "re " << re << endl;
	outStream << "nu " << nu << endl;
	outStream << "minRho " << minRho << endl;
	outStream << "rhoS " << rhoS << endl;
	
	Simulation_Fluid::_outputSettings(outStream);
}

void Sim_Multiphase::_inputSettings(istream& inStream)
{
	string variableName;
	
	inStream >> variableName;
	if (variableName != "Multiphase")
	{
		cout << "Error in deserialization - Simulation_Multiphase\n";
		abort();
	}
	
	// read data
	inStream >> variableName;
	assert(variableName=="re");
	inStream >> re;
	inStream >> variableName;
	assert(variableName=="nu");
	inStream >> nu;
	inStream >> variableName;
	assert(variableName=="minRho");
	inStream >> minRho;
	inStream >> variableName;
	assert(variableName=="rhoS");
	inStream >> rhoS;
	
	Simulation_Fluid::_inputSettings(inStream);
}

Sim_Multiphase::Sim_Multiphase(const int argc, const char ** argv) : Simulation_Fluid(argc, argv), gravity{0,-9.81}, dtCFL(0), dtFourier(0), re(0), nu(0), minRho(0), rhoS(1), bSplit(false)
{
#ifdef _MULTIGRID_
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
#endif // _MULTIGRID_
	
	if (rank==0)
	{
		cout << "====================================================================================================================\n";
		cout << "\t\t\tMultiphase flow\n";
		cout << "====================================================================================================================\n";
	}
}

Sim_Multiphase::~Sim_Multiphase()
{
}

void Sim_Multiphase::init()
{
	Simulation_Fluid::init();
	
	if (!bRestart)
	{
		// simulation settings
		bSplit = parser("-split").asBool(false);
		nu = parser("-nu").asDouble(1e-2);
		rhoS = parser("-rhoS").asDouble(1);
		minRho = min((Real)1.,(Real)rhoS);
		
		stringstream ss;
		ss << path2file << "_settings.dat";
		ofstream myfile(ss.str(), fstream::app);
		_dumpSettings(cout);
		_dumpSettings(myfile);
		
		if (rank==0)
			if (bSplit)
				cout << "Using split method with constant coefficients Poisson solver\n";
			else
				cout << "Solving full variable coefficient Poisson equation for pressure\n";
		
		_ic();
	}
	
	pipeline.clear();
	pipeline.push_back(new CoordinatorCleanTmp(grid));
	pipeline.push_back(new CoordinatorGravity(gravity, grid));
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid));
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
	pipeline.push_back(new CoordinatorUpdate(grid));
	pipeline.push_back(new CoordinatorPressure<Lab>(minRho, &step, bSplit, grid, rank, nprocs));
	
	if (rank==0)
	{
		cout << "Coordinator/Operator ordering:\n";
		for (int c=0; c<pipeline.size(); c++)
			cout << "\t" << pipeline[c]->getName() << endl;
	}
}

void Sim_Multiphase::simulate()
{
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
	
	double vOld = 0;
	
#ifdef _MULTIGRID_
	MPI_Barrier(MPI_COMM_WORLD);
#endif // _MULTIGRID_
	double nextDumpTime = time;
	double maxU = 0;
	
	while (true)
	{
		if (rank==0)
		{
			vector<BlockInfo> vInfo = grid->getBlocksInfo();
			
			// choose dt (CFL, Fourier)
			profiler.push_start("DT");
			maxU = findMaxUOMP(vInfo,*grid);
			dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
			dtCFL     = maxU==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/abs(maxU);
			assert(!std::isnan(maxU));
			dt = min(dtCFL,dtFourier);
			if (dumpTime>0)
				dt = min(dt,nextDumpTime-_nonDimensionalTime());
			if (endTime>0)
				dt = min(dt,endTime-_nonDimensionalTime());
			if (verbose)
				cout << "dt (Fourier, CFL): " << dt << " " << dtFourier << " " << dtCFL << endl;
			profiler.pop_stop();
		}
		MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		if (dt!=0)
		{
			for (int c=0; c<pipeline.size(); c++)
			{
#ifdef _MULTIGRID_
				MPI_Barrier(MPI_COMM_WORLD);
#endif // _MULTIGRID_
				profiler.push_start(pipeline[c]->getName());
				if (rank == 0 || pipeline[c]->getName()=="Pressure")
					(*pipeline[c])(dt);
				profiler.pop_stop();
			}
			
			time += dt;
			step++;
		}
		
		if (rank==0)
		{
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
			
			if (step % 100 == 0)
				profiler.printSummary();
		}
		
		// check nondimensional time
		if ((endTime>0 && abs(_nonDimensionalTime()-endTime) < 10*std::numeric_limits<Real>::epsilon()) || (nsteps!=0 && step>=nsteps))
		{
			if (rank==0)
			{
				profiler.push_start("Dump");
				stringstream ss;
				ss << path2file << "-Final.vti";
				cout << ss.str() << endl;
				
				dumper.Write(*grid, ss.str());
				
				vector<BlockInfo> vInfo = grid->getBlocksInfo();
				Layer vorticity(sizeX,sizeY,1);
				processOMP<Lab, OperatorVorticity>(vorticity,vInfo,*grid);
				stringstream sVort;
				sVort << path2file << "Vorticity-Final.vti";
				dumpLayer2VTK(step,sVort.str(),vorticity,1);
				profiler.pop_stop();
				
				profiler.printSummary();
			}
			exit(0);
		}
	}
}