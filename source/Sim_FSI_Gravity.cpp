//
//  Sim_FSI_Gravity.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/26/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Sim_FSI_Gravity.h"

#include "ProcessOperatorsOMP.h"

#include "CoordinatorIC.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPenalization.h"
#include "CoordinatorComputeShape.h"
#include "CoordinatorPressure.h"
#include "CoordinatorGravity.h"
#include "CoordinatorBodyVelocities.h"
#include "CoordinatorTankTread.h"

void Sim_FSI_Gravity::_diagnostics()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	double drag = 0;
	double volS = 0;
	double volF = 0;
	double pMin = 10;
	double pMax = 0;
	double centerMassX = 0;
	double centerMassY = 0;
	double centroidX = 0;
	double centroidY = 0;
	double mass = 0;
	double volume = 0;
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for schedule(static) reduction(+:drag) reduction(+:volS) reduction(+:volF) reduction(max:pMax) reduction (min:pMin) reduction(+:centerMassX) reduction(+:centerMassY) reduction(+:centroidX) reduction(+:centroidY) reduction(+:mass) reduction(+:volume)
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
		{
			double p[2] = {0,0};
			info.pos(p, ix, iy);
			const double chi = b(ix,iy).chi;
			const double rhochi = b(ix,iy).rho * chi;
			centerMassX += p[0] * rhochi;
			centerMassY += p[1] * rhochi;
			centroidX += p[0] * chi;
			centroidY += p[1] * chi;
			mass += rhochi;
			volume += chi;
			
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
	
	centerMassX /= mass;
	centerMassY /= mass;
	centroidX /= volume;
	centroidY /= volume;
	const double rhoSAvg = mass/volume;
	
	double cD = 2*drag/(uBody[1]*uBody[1]*shape->getCharLength());
	cD = abs(uBody[1])>0 ? cD : 1e10;
	const Real Re_uBody = shape->getCharLength()*sqrt(uBody[0]*uBody[0]+uBody[1]*uBody[1])/nu;
	Real centroid[2], com[2];
	shape->getCentroid(centroid);
	shape->getCenterOfMass(com);
    
    Real dragPD = (dragP[1]+dragV)*dh*dh;
    Real cD_PD = 2*dragPD/(uBody[1]*uBody[1]*shape->getCharLength());
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	if (verbose)
		cout << step << " " << time << " " << dt << " " << bpdx << " " << lambda << " " << cD << " " << Re_uBody << " " << centroid[0] << " " << centroid[1] << " " << uBody[0] << " " << uBody[1] << " " << shape->getOrientation() << " " << omegaBody << " " << com[0] << " " << com[1] << " " << centerMassX << " " << centerMassY << " " << centroidX << " " << centroidY << " " << rhoSAvg << " " << cD_PD << endl;
	myfile << step << " " << time << " " << dt << " " << bpdx << " " << lambda << " " << cD << " " << Re_uBody << " " << centroid[0] << " " << centroid[1] << " " << uBody[0] << " " << uBody[1] << " " << shape->getOrientation() << " " << omegaBody << " " << com[0] << " " << com[1] << " " << centerMassX << " " << centerMassY << " " << centroidX << " " << centroidY << " " << rhoSAvg << " " << cD_PD << endl;
}

void Sim_FSI_Gravity::_dumpSettings(ostream& outStream)
{
#ifdef _MULTIGRID_
	if (rank==0)
#endif // _MULTIGRID_
	{
		outStream << "--------------------------------------------------------------------\n";
		outStream << "Physical Settings\n";
		outStream << "\tradius\t" << shape->getCharLength()*.5 << endl;
		outStream << "\tnu\t" << nu << endl;
		outStream << "\tMinRhoS\t" << shape->getMinRhoS() << endl;
		Real center[2];
		shape->getCentroid(center);
		outStream << "\tyPos\t" << center[1] << endl;
		
		outStream << "\nSimulation Settings\n";
		outStream << "\tnsteps\t" << nsteps << endl;
		outStream << "\tTend\t" << endTime << endl;
		outStream << "\tlambda\t" << lambda << endl;
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

void Sim_FSI_Gravity::_outputSettings(ostream &outStream)
{
	outStream << "Gravity_FSI\n";
	outStream << "uBody " << uBody[0] << endl;
	outStream << "vBody " << uBody[1] << endl;
	outStream << "omegaBody " << omegaBody << endl;
	outStream << "re " << re << endl;
	outStream << "nu " << nu << endl;
	
	Simulation_FSI::_outputSettings(outStream);
}

void Sim_FSI_Gravity::_inputSettings(istream& inStream)
{
	string variableName;
	
	inStream >> variableName;
	if (variableName != "Gravity_FSI")
	{
		cout << "Error in deserialization - Simulation_Gravity_FSI\n";
		abort();
	}
	
	// read data
	inStream >> variableName;
	assert(variableName=="uBody");
	inStream >> uBody[0];
	inStream >> variableName;
	assert(variableName=="vBody");
	inStream >> uBody[1];
	inStream >> variableName;
	assert(variableName=="omegaBody");
	inStream >> omegaBody;
	inStream >> variableName;
	assert(variableName=="re");
	inStream >> re;
	inStream >> variableName;
	assert(variableName=="nu");
	inStream >> nu;
	
	Simulation_FSI::_inputSettings(inStream);
}

Sim_FSI_Gravity::Sim_FSI_Gravity(const int argc, const char ** argv) : Simulation_FSI(argc, argv), uBody{0,0}, omegaBody(0), gravity{0,-9.81}, dtCFL(0), dtFourier(0), dtBody(0), re(0), nu(0), minRho(0), bSplit(false)
{
#ifdef _MULTIGRID_
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
	if (rank!=0)
		omp_set_num_threads(1);
#endif // _MULTIGRID_
	
	if (rank==0)
	{
		cout << "====================================================================================================================\n";
		cout << "\t\t\tFlow past a falling cylinder\n";
		cout << "====================================================================================================================\n";
	}
}

Sim_FSI_Gravity::~Sim_FSI_Gravity()
{
}

void Sim_FSI_Gravity::init()
{
	Simulation_FSI::init();
	
	if (!bRestart)
	{
		// simulation settings
		bSplit = parser("-split").asBool(false);
		nu = parser("-nu").asDouble(1e-2);
		minRho = min((Real)1.,shape->getMinRhoS());
		
		gravity[1] = -parser("-g").asDouble(9.81);
		
		const Real aspectRatio = (Real)(bpdx*FluidBlock::sizeX)/(Real)(bpdy*FluidBlock::sizeY);
		Real center[2] = {parser("-xpos").asDouble(.5*aspectRatio),parser("-ypos").asDouble(.85)};
		shape->setCentroid(center);
		
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
#ifndef _MULTIPHASE_
	pipeline.push_back(new CoordinatorAdvection<Lab>(&uBody[0], &uBody[1], grid));
#else
	pipeline.push_back(new CoordinatorAdvection<Lab>(&uBody[0], &uBody[1], grid, 1));
#endif
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, &uBody[0], &uBody[1], &dragV, grid));
	pipeline.push_back(new CoordinatorGravity(gravity, grid));
#ifdef _TANKTREADING_
	pipeline.push_back(new CoordinatorTankTread(shape,grid));
#endif
	pipeline.push_back(new CoordinatorPressure<Lab>(minRho, gravity, &uBody[0], &uBody[1], &dragP[0], &dragP[1], &step, bSplit, grid, rank, nprocs));
	pipeline.push_back(new CoordinatorBodyVelocities(&uBody[0], &uBody[1], &omegaBody, shape, &lambda, grid));
	pipeline.push_back(new CoordinatorPenalization(&uBody[0], &uBody[1], &omegaBody, shape, &lambda, grid));
#ifndef _MOVING_FRAME_
	pipeline.push_back(new CoordinatorComputeShape(&uBody[0], &uBody[1], &omegaBody, shape, grid));
#endif
	
	if (rank==0)
	{
		cout << "Coordinator/Operator ordering:\n";
		for (int c=0; c<pipeline.size(); c++)
			cout << "\t" << pipeline[c]->getName() << endl;
	}
}

void Sim_FSI_Gravity::simulate()
{
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
	
	double vOld = 0;
	
#ifdef _MULTIGRID_
	MPI_Barrier(MPI_COMM_WORLD);
#endif // _MULTIGRID_
	double nextDumpTime = time;
	double maxU = 0;
	double maxA = 0;
	
	while (true)
	{
		if (rank==0)
		{
			vector<BlockInfo> vInfo = grid->getBlocksInfo();
			
			// choose dt (CFL, Fourier)
			profiler.push_start("DT");
#ifndef _MOVING_FRAME_
			maxU = findMaxUOMP(vInfo,*grid);
#else
			maxU = findMaxUOMP(uBody[0],uBody[1],vInfo,*grid);
#endif
#ifdef _MULTIPHASE_
			dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu*min(shape->getMinRhoS(),(Real)1);
#else
#ifdef _CONSTNU_
			dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
#else
			dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu*min(shape->getMinRhoS(),(Real)1);
#endif
#endif
			dtCFL     = maxU==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/abs(maxU);
			dtBody    = max(abs(uBody[0]),abs(uBody[1]))==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/max(abs(uBody[0]),abs(uBody[1]));
			assert(!std::isnan(maxU));
			assert(!std::isnan(maxA));
			assert(!std::isnan(uBody[0]));
			assert(!std::isnan(uBody[1]));
			dt = min(min(dtCFL,dtFourier),dtBody);
#ifdef _PARTICLES_
			maxA = findMaxAOMP<Lab>(vInfo,*grid);
			dtLCFL = maxA==0 ? 1e5 : LCFL/abs(maxA);
			dt = min(dt,dtLCFL);
#endif
			if (dumpTime>0)
				dt = min(dt,nextDumpTime-_nonDimensionalTime());
			if (endTime>0)
				dt = min(dt,endTime-_nonDimensionalTime());
			if (verbose)
				cout << "time, dt (Fourier, CFL, body): " << time << " " << dt << " " << dtFourier << " " << dtCFL << " " << dtBody << endl;
			profiler.pop_stop();
		}
#ifdef _MULTIGRID_
		MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif // _MULTIGRID_
		
		if (dt!=0)
		{
#ifdef _DLM_
			lambda = dlm/dt;
#endif
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
			//if (step<100)
			{
				// this still needs to be corrected to the frame of reference!
				// this test only works for constant density disks as it is written now
				double accM = (uBody[1]-vOld)/dt;
				vOld = uBody[1];
				double accT = (shape->getMinRhoS()-1)/(shape->getMinRhoS()+1) * gravity[1];
				double accN = (shape->getMinRhoS()-1)/(shape->getMinRhoS()  ) * gravity[1];
				if (verbose) cout << "Acceleration with added mass (measured, expected, no added mass)\t" << accM << "\t" << accT << "\t" << accN << endl;
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
				profiler.pop_stop();
				
				profiler.printSummary();
			}
			exit(0);
		}
	}
}