//
//  Sim_FSI_Gravity.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/26/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Sim_FSI_Gravity.h"

#include "ProcessOperatorsOMP.h"
//#include "OperatorIC.h"
//#include "OperatorAdvection.h"
//#include "OperatorDiffusion.h"
//#include "OperatorPenalization.h"
//#include "OperatorDivergence.h"
#include "OperatorVorticity.h"
//#include "PoissonSolverScalarFFTW.h"
//#include "OperatorGradP.h"
//#include "OperatorComputeShape.h"
//#include "OperatorGravity.h"
//#include "OperatorSplitP.h"

#include "CoordinatorIC.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPenalization.h"
#include "CoordinatorComputeShape.h"
#include "CoordinatorPressure.h"
#include "CoordinatorGravity.h"
#include "CoordinatorBodyVelocities.h"

void Sim_FSI_Gravity::_dumpDivergence(const int step, const Real rho0, const Real dt)
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const int size = bpdx * FluidBlock::sizeX;
	Layer divergence(size,size,1);
	Layer divergenceSplit(size,size,1);
	processOMP<Lab, OperatorDivergenceLayer>(divergence, vInfo, *grid);
	processOMP<Lab, OperatorDivergenceSplitLayer>(divergenceSplit, rho0, dt, step, vInfo, *grid);
	
	stringstream sDiv, sDivSplit;
	sDiv << path2file << "Divergence-" << step << ".vti";
	sDivSplit << path2file << "DivergenceSplit-" << step << ".vti";
	dumpLayer2VTK(step,sDiv.str(),divergence,1);
	dumpLayer2VTK(step,sDivSplit.str(),divergenceSplit,1);
}

void Sim_FSI_Gravity::_diagnostics()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	double drag = 0;
	double volS = 0;
	double volF = 0;
	double pMin = 10;
	double pMax = 0;
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for schedule(static) reduction(+:drag) reduction(+:volS) reduction(+:volF) reduction(max:pMax) reduction (min:pMin)
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				pMin = min(pMin,(double)b(ix,iy).p);
				pMax = max(pMax,(double)b(ix,iy).p);
				
				if (b(ix,iy).chi>0)
				{
					drag += (b(ix,iy).v-uBody[1]) * b(ix,iy).chi; // this depends on the direction of movement - here vertical!
					volS += b(ix,iy).chi;
					volF += (1-b(ix,iy).chi);
				}
				
				if (std::isnan(b(ix,iy).u) ||
					std::isnan(b(ix,iy).v) ||
					std::isnan(b(ix,iy).rho) ||
					std::isnan(b(ix,iy).chi) ||
					std::isnan(b(ix,iy).p))
				{
					cout << "NaN Error - Aborting now!\n";
					abort();
				}
			}
	}
	
	drag *= dh*dh*lambda;
	volS *= dh*dh;
	volF *= dh*dh;
	
	double cD = 2*drag/(uBody[1]*uBody[1]*shape->getCharLength());
	cD = abs(uBody[1])>0 ? cD : 1e10;
	const Real Re_uBody = shape->getCharLength()*abs(uBody[1])/nu;
	Real center[2];
	shape->getPosition(center);
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	if (verbose)
		cout << step << " " << time << " " << dt << " " << bpdx << " " << lambda << " " << cD << " " << Re_uBody << " " << center[0] << " " << center[1] << " " << uBody[0] << " " << uBody[1] << " " << shape->getOrientation() << endl;
	myfile << step << " " << time << " " << dt << " " << bpdx << " " << lambda << " " << cD << " " << Re_uBody << " " << center[0] << " " << center[1] << " " << uBody[0] << " " << uBody[1] << " " << shape->getOrientation() << endl;
}

void Sim_FSI_Gravity::_dumpSettings(ostream& mystream)
{
#ifdef _MULTIGRID_
	if (rank==0)
#endif // _MULTIGRID_
	{
		mystream << "--------------------------------------------------------------------\n";
		mystream << "Physical Settings\n";
		mystream << "\tradius\t" << shape->getCharLength()*.5 << endl;
		mystream << "\tnu\t" << nu << endl;
		mystream << "\tRhoS\t" << shape->getRhoS() << endl;
		Real center[2];
		shape->getPosition(center);
		mystream << "\tyPos\t" << center[1] << endl;
		
		mystream << "\nSimulation Settings\n";
		mystream << "\tnsteps\t" << nsteps << endl;
		mystream << "\tTend\t" << endTime << endl;
		mystream << "\tlambda\t" << lambda << endl;
#ifdef _MULTIGRID_
		mystream << "\tsplit\t" << (bSplit ? "true" : "false") << endl;
#endif // _MULTIGRID_
		mystream << "\tCFL\t" << CFL << endl;
		mystream << "\tpath2file\t" << path2file << endl;
		mystream << "\tbpdx\t" << bpdx << endl;
		mystream << "\tbpdy\t" << bpdy << endl;
#ifdef _PERIODIC_
		mystream << "\tBC\t\tperiodic\n";
#else // _PERIODIC_
		mystream << "\tBC\tmixed\n";
#endif // _PERIODIC_
#ifdef _MULTIGRID_
		mystream << "\tPoisson\tMultigrid\n";
#endif // _MULTIGRID_
#ifdef _SPLIT_
		mystream << "\tPoisson\tFFTW Split\n";
#endif // _SPLIT_
		mystream << "--------------------------------------------------------------------\n";
	}
}

void Sim_FSI_Gravity::_ic()
{
#ifdef _MULTIGRID_
	if (rank==0)
#endif // _MULTIGRID_
	{
		// setup initial conditions
		CoordinatorIC coordIC(shape,0,grid);
		profiler.push_start(coordIC.getName());
		coordIC(0);
		
		stringstream ss;
		ss << path2file << "-IC.vti";
		dumper.Write(*grid, ss.str());
		profiler.pop_stop();
	}
}

double Sim_FSI_Gravity::_nonDimensionalTime()
{
	return time; // how to nondimensionalize here? based on Galileo number?
}

Sim_FSI_Gravity::Sim_FSI_Gravity(const int argc, const char ** argv) : Simulation_FSI(argc, argv), uBody{0,0}, omegaBody(0), gravity{0,-9.81}, dtCFL(0), dtFourier(0), dtBody(0), re(0), nu(0), CFL(0), minRho(0), bSplit(false)
#ifdef _MULTIGRID_
, rank(0), nprocs(0)
#endif // _MULTIGRID_
{
#ifdef _MULTIGRID_
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
	if (rank==0)
#endif // _MULTIGRID_
	{
		cout << "====================================================================================================================\n";
		cout << "\t\t\tFlow past a falling cylinder\n";
		cout << "====================================================================================================================\n";
	}
	
	// simulation settings
	bSplit = parser("-split").asBool(false);
	nu = parser("-nu").asDouble(1e-2);
	CFL = parser("-CFL").asDouble(.1);
	minRho = min((Real)1.,shape->getRhoS());
	
	const float aspectRatio = (float)bpdx/(float)bpdy;
	Real center[2] = {.5*aspectRatio,parser("-ypos").asDouble(.85)};
	shape->setPosition(center);
	
	stringstream ss;
	ss << path2file << "_settings.dat";
	ofstream myfile(ss.str(), fstream::app);
	_dumpSettings(cout);
	_dumpSettings(myfile);
	
#ifdef _MULTIGRID_
	if (rank==0)
		if (bSplit)
			cout << "Using split method with constant coefficients Poisson solver\n";
		else
			cout << "Solving full variable coefficient Poisson equation for pressure\n";
#endif // _MULTIGRID_
	
	_ic();
	
	pipeline.clear();
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid));
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
	pipeline.push_back(new CoordinatorPenalization(&uBody[0], &uBody[1], &omegaBody, shape, lambda, grid));
	pipeline.push_back(new CoordinatorComputeShape(&uBody[0], &uBody[1], &omegaBody, shape, grid));
#ifdef _MULTIGRID_
	pipeline.push_back(new CoordinatorPressure<Lab>(rank, nprocs, minRho, &step, bSplit, grid));
#else // _MULTIGRID_
	pipeline.push_back(new CoordinatorPressure<Lab>(0, 1, minRho, &step, bSplit, grid));
#endif // _MULTIGRID_
	pipeline.push_back(new CoordinatorGravity(gravity, grid));
	pipeline.push_back(new CoordinatorBodyVelocities(&uBody[0], &uBody[1], &omegaBody, lambda, grid));
	
#ifdef _MULTIGRID_
	if (rank==0)
#endif // _MULTIGRID_
	{
	cout << "Coordinator/Operator ordering:\n";
	for (int c=0; c<pipeline.size(); c++)
		cout << "\t" << pipeline[c]->getName() << endl;
	}
	
#ifdef _MULTIGRID_
	MPI_Barrier(MPI_COMM_WORLD);
#endif // _MULTIGRID_
}

Sim_FSI_Gravity::~Sim_FSI_Gravity()
{
}

void Sim_FSI_Gravity::simulate()
{
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
	
	double vOld = 0;
	Real oldAccVort[2] = {0,0};
	
#ifdef _MULTIGRID_
	MPI_Barrier(MPI_COMM_WORLD);
#endif // _MULTIGRID_
	time = 0;
	double nextDumpTime = 0;
	double maxU = 0;
	while (true)
	{
#ifdef _MULTIGRID_
		if (rank==0)
#endif // _MULTIGRID_
		{
			//cout << "step " << step << endl;
			vector<BlockInfo> vInfo = grid->getBlocksInfo();
			
			// choose dt (CFL, Fourier)
			profiler.push_start("DT");
			maxU = findMaxUOMP(vInfo,*grid);
			dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
			dtCFL     = maxU==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/abs(maxU);
			dtBody    = max(abs(uBody[0]),abs(uBody[1]))==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/max(abs(uBody[0]),abs(uBody[1]));
			dt = min(min(dtCFL,dtFourier),dtBody);
			if (dumpTime>0)
				dt = min(dt,nextDumpTime-_nonDimensionalTime());
			if (endTime>0)
				dt = min(dt,endTime);
			if (verbose)
				cout << "dt (Fourier, CFL, body): " << dt << " " << dtFourier << " " << dtCFL << " " << dtBody << endl;
			profiler.pop_stop();
			
			// gravity
			profiler.push_start(pipeline[5]->getName());
			(*pipeline[5])(dt);
			profiler.pop_stop();
			
			// pressure
			profiler.push_start(pipeline[4]->getName());
		}
		(*pipeline[4])(dt);
#ifdef _MULTIGRID_
		if (rank==0)
#endif // _MULTIGRID_
		{
			profiler.pop_stop();
			
			vector<BlockInfo> vInfo = grid->getBlocksInfo();
			
			// advection
			profiler.push_start(pipeline[0]->getName());
			(*pipeline[0])(dt);
			profiler.pop_stop();
			
			// diffusion
			profiler.push_start(pipeline[1]->getName());
			(*pipeline[1])(dt);
			profiler.pop_stop();
			
			if (step>50)
			{
				profiler.push_start("Body");
				vOld = uBody[1];
				(*pipeline[6])(dt);
				//computeBodyVelocity(vInfo, *grid, uBody, omegaBody, shape->getRhoS(), gravity, dt, lambda);
				(*pipeline[3])(dt);
				profiler.pop_stop();
			}
			
			// penalization
			profiler.push_start(pipeline[2]->getName());
			(*pipeline[2])(dt);
			profiler.pop_stop();
			
			time += dt;
			step++;
			
			if (step<100)
			{
				// this still needs to be corrected to the frame of reference!
				double accM = (uBody[1]-vOld)/dt;
				double accT = (shape->getRhoS()-1)/(shape->getRhoS()+1) * gravity[1];
				double accN = (shape->getRhoS()-1)/(shape->getRhoS()) * gravity[1];
				cout << "Acceleration with added mass (measured, expected, no added mass)\t" << accM << "\t" << accT << "\t" << accN << endl;
				stringstream ss;
				ss << path2file << "_addedmass.dat";
				ofstream myfile(ss.str(), fstream::app);
				myfile << step << " " << time << " " << accM << " " << accT << " " << accN << endl;
			}
			
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
}