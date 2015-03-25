//
//  Sim_FSI_Gravity.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/26/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Sim_FSI_Gravity.h"
#include "ProcessOperatorsOMP.h"
#include "InterfaceFortran.h"
#include "Timer.h"
#include "LayerToVTK.h"
#include <sstream>
#include <cmath>

void Sim_FSI_Gravity::_solvePressure()
{
	// need an interface that is the same for all solvers - this way the defines can be removed more cleanly
#ifdef _MULTIGRID_
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	// pressure
#ifdef _SPLIT_
#ifdef _SP_COMP_
	PoissonSolverScalarFFTW<FluidGrid, StreamerDiv> pressureSolver(NTHREADS);
	processOMP<Lab, OperatorDivergenceSplit>(dt, minRho, step, vInfo, *grid);
	pressureSolver.solve(*grid,false);
	processOMP<Lab, OperatorGradPSplit>(dt, minRho, step, vInfo, *grid);
#else
	cout << "FFTW double precision not supported - aborting now!\n";
	abort();
#endif // _SP_COMP_
#endif // _SPLIT_
#ifdef _FORTRAN_SOLVER_
	processOMP<Lab, OperatorDivergenceSplit>(dt, minRho, step, vInfo, *grid);
	
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
	float startx = 0;
	float endx = 1;
	int m = sizeX-1;
	int bcx = 0;
	float starty = 0;
	float endy = 1;
	int n = sizeY-1;
	int bcy = 4; // dirichlet at y=1 and Neumann at y=0
	float * bdc = new float[sizeY];
	for (int i=0; i<sizeY; i++) bdc[i] = 0;
	float lambda = 0;
	float * grhs = new float[sizeX*sizeY];
	// fill arrays
	int idmn = sizeX;
	float pertrb = 0;
	int ierror = 0;
	
	processOMP_preparePoissonFortran(grhs,sizeY,vInfo,*grid);
	
	hwscrt_(&startx, &endx, &m, &bcx, NULL, NULL,
			&starty, &endy, &n, &bcy, bdc, NULL,
			&lambda, grhs, &idmn, &pertrb, &ierror);
	
	if (ierror!=0)
		cout << "ierror " << ierror << endl;
	
	processOMP_readPoissonFortran(grhs,pertrb,ierror,sizeY,vInfo,*grid);
	delete [] bdc;
	delete [] grhs;
	
	processOMP<Lab, OperatorGradPSplit>(dt, minRho, step, vInfo, *grid);
#endif
#ifdef _JACOBI_
	processOMP<Lab, OperatorDivergence>(dt, vInfo, *grid);
	processOMP_Jacobi<Lab, OperatorVarCoeffPoisson>(vInfo,*grid);
	processOMP<Lab, OperatorGradP>(dt, vInfo, *grid);
#endif
#ifdef _MULTIGRID_
	if (rank==0)
		if (bSplit)
			processOMP<Lab, OperatorDivergenceSplit>(dt, minRho, step, vInfo, *grid);
		else
			processOMP<Lab, OperatorDivergence>(dt, vInfo, *grid);
	mg.setup(grid, bSplit, rank, nprocs);
	mg();
	if (rank==0)
		if (bSplit)
			processOMP<Lab, OperatorGradPSplit>(dt, minRho, step, vInfo, *grid);
		else
			processOMP<Lab, OperatorGradP>(dt, vInfo, *grid);
#endif
	
#ifdef _MULTIGRID_
	if (rank==0)
#endif
		updatePressuresOMP(vInfo, *grid);
}

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
	myfile << step << " " << time << " " << dt << " " << bpdx << " " << lambda << " " << cD << " " << Re_uBody << " " << center[1] << " " << uBody[1] << " " << shape->getOrientation() << endl;
	// how to nondimensionalize time for an accelerating body?
}

void Sim_FSI_Gravity::_dumpSettings(ostream& mystream)
{
#ifdef _MULTIGRID_
	if (rank==0)
#endif
	{
		mystream << "--------------------------------------------------------------------\n";
		mystream << "Physical Settings\n";
		mystream << "\tradius\t" << shape->getCharLength()*.5 << endl;
		mystream << "\tnu\t" << nu << endl;
		mystream << "\tRhoS\t" << rhoS << endl;
		Real center[2];
		shape->getPosition(center);
		mystream << "\tyPos\t" << center[1] << endl;
		
		mystream << "\nSimulation Settings\n";
		mystream << "\tnsteps\t" << nsteps << endl;
		mystream << "\tTend\t" << tEnd << endl;
		mystream << "\tlambda\t" << lambda << endl;
#ifdef _MULTIGRID_
		mystream << "\tsplit\t" << (bSplit ? "true" : "false") << endl;
#endif
		mystream << "\tCFL\t" << CFL << endl;
		mystream << "\tpath2file\t" << path2file << endl;
		mystream << "\tbpdx\t" << bpdx << endl;
		mystream << "\tbpdy\t" << bpdy << endl;
#ifdef _PERIODIC_
		mystream << "\tBC\t\tperiodic\n";
#else
		mystream << "\tBC\tmixed\n";
#endif
#ifdef _MULTIGRID_
		mystream << "\tPoisson\tMultigrid\n";
#endif
#ifdef _SPLIT_
		mystream << "\tPoisson\tFFTW Split\n";
#endif
#ifdef _JACOBI_
		mystream << "\tPoisson\tJacobi\n";
#endif
#ifdef _FORTRAN_SOLVER_
		mystream << "\tPoisson\tFishPack\n";
#endif
		mystream << "--------------------------------------------------------------------\n";
	}
}

void Sim_FSI_Gravity::setup()
{
	// simulation settings
	nsteps = parser("-nsteps").asInt(1000000);
	uinf = 0;//parser("-uinf").asDouble(0.1);
	re = 0;//parser("-Re").asDouble(100);
	nu = parser("-nu").asDouble(1e-2);
	lambda = parser("-lambda").asDouble(1e8);
	tEnd = parser("-tend").asDouble(50);
	rhoS = parser("-rhoS").asDouble(1.1);
	CFL = parser("-CFL").asDouble(.1);
	minRho = min(1.,rhoS);
	
	// output settings
	path2file = parser("-file").asString("../data/FlowPastFallingCylinder");
	
	// computational space settings
	parser.set_strict_mode();
	bpdx = parser("-bpdx").asInt();
	bpdy = parser("-bpdy").asInt();
	parser.unset_strict_mode();
	
	const float aspectRatio = (float)bpdx/(float)bpdy;
	Real centerOfMass[2] = {.5*aspectRatio,parser("-ypos").asDouble(.85)};
	bool bPeriodic[2] = {false,false};
	string shapeType = parser("-shape").asString("disk");
	
	if (shapeType=="disk")
	{
		Real radius = parser("-radius").asDouble(0.1);
		shape = new Disk(centerOfMass, radius, rhoS, 2, 2, bPeriodic);
	}
	else if (shapeType=="ellipse")
	{
		Real semiAxis[2] = {parser("-semiAxisX").asDouble(0.1),parser("-semiAxisY").asDouble(0.2)};
		Real angle = parser("-angle").asDouble(0.0);
		shape = new Ellipse(centerOfMass, semiAxis, angle, rhoS, 2, 2, bPeriodic);
	}
		
	
	stringstream ss;
	ss << path2file << "_settings.dat";
	ofstream myfile(ss.str(), fstream::app);
	_dumpSettings(cout);
	_dumpSettings(myfile);
	
#ifdef _PERIODIC_
	if (bpdx!=bpdy)
	{
		cout << "Assuming square domain for periodicity in OperatorComputeDisk\n";
		abort();
	}
#endif
	
#ifdef _MULTIGRID_
	if (rank==0)
		if (bSplit)
			cout << "Using split method with constant coefficients Poisson solver\n";
		else
			cout << "Solving full variable coefficient Poisson equation for pressure\n";
#endif
	
	grid = new FluidGrid(bpdx,bpdy,1);
	assert(grid != NULL);
	
	Timer timerIC;
#ifdef _MULTIGRID_
	if (rank==0)
#endif
	{
		// setup initial conditions
		timerIC.start();
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		OperatorIC ic(shape, 0);
		FluidBlockProcessing::process(vInfo, ic, grid);
		
		dt = 1e-8;
		
		// this is required to create the initial pressure gradient
		//processOMP_hydrostaticTerm<OperatorGravity>(gravity, dt, rhoS, vInfo, *grid);
		//processOMP<Lab, OperatorGradP>(dt, vInfo, *grid);
		
	}
	//_solvePressure(); //  put mgh instead of solving!
#ifdef _MULTIGRID_
	if (rank==0)
#endif
	{
		stringstream ss;
		ss << path2file << "-IC.vti" ;
		cout << ss.str() << endl;
		
		dumper.Write(*grid, ss.str());
		
		//_dumpDivergence(-1,minRho,dt);
		//_diagnostics();
		
		double timeIC = timerIC.stop();
		cout << "Time IC (u, p):\t" << timeIC << endl;// << "\t" << timePoissonIC << endl;
	}
	
#ifdef _MULTIGRID_
	MPI_Barrier(MPI_COMM_WORLD);
	//MPI_Abort(MPI_COMM_WORLD,1);
#endif
}

void Sim_FSI_Gravity::simulate()
{
	time = 0;
	double nextDumpTime = 0;
	double maxU = uinf;
	Timer timer, timerPressure;
	
	double timeDT = 0;
	double timeAdvection = 0;
	double timeDiffusion = 0;
	double timePressure = 0;
	double timeGravity = 0;
	double timeUB = 0;
	double timePenalization = 0;
	double timeDiagnostics = 0;
	double timeDump = 0;
	
	double timeDivergence = 0;
	double timePoisson = 0;
	double timeUpdatePressure = 0;
	
	double vOld = 0;
	Real oldAccVort[2] = {0,0};
	
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
	
#ifdef _MULTIGRID_
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	for(step=0; step<nsteps; ++step)
	{
#ifdef _MULTIGRID_
		if (rank==0)
#endif
		{
			//cout << "step " << step << endl;
			vector<BlockInfo> vInfo = grid->getBlocksInfo();
			
			// choose dt (CFL, Fourier)
			timer.start();
			maxU = findMaxUOMP(vInfo,*grid);
			dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
			dtCFL     = maxU==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/abs(maxU);
			dtBody    = max(abs(uBody[0]),abs(uBody[1]))==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/max(abs(uBody[0]),abs(uBody[1]));
			dt = min(min(dtCFL,dtFourier),dtBody);
			//cout << "dt (Fourier, CFL, body): " << dt << " " << dtFourier << " " << dtCFL << " " << dtBody << endl;
			timeDT += timer.stop();
			
			//lambda = 1./dt;
			
			// gravity
			timer.start();
			processOMP_hydrostaticTerm<OperatorGravity>(gravity, dt, rhoS, vInfo, *grid);
			timeGravity += timer.stop();
			
			// pressure
			timer.start();
		}
		//timerPressure.start();
		_solvePressure();
		//timeUpdatePressure += timerPressure.stop();
#ifdef _MULTIGRID_
		if (rank==0)
#endif
		{
			timePressure += timer.stop();
			
			vector<BlockInfo> vInfo = grid->getBlocksInfo();
			
			// advection
			timer.start();
			resetOMP(vInfo, *grid);
			processOMP< Lab,OperatorAdvection<Mp4> >(dt,vInfo,*grid);
			updateOMP(vInfo, *grid);
			timeAdvection += timer.stop();
			
			// diffusion
			timer.start();
			resetOMP(vInfo, *grid);
			processOMP<Lab,OperatorDiffusion>(dt,nu,vInfo,*grid);
			updateOMP(vInfo, *grid);
			timeDiffusion += timer.stop();
			
			//_dumpDivergence(step,minRho,dt);
			
			// this needs to be debugged, maybe first go back to the original version, instead of the optimized one
			//Layer vortTmp(sizeX,sizeY,1);
			//processOMP<Lab, OperatorVorticity>(vortTmp,vInfo,*grid);
			//computeForcesFromVorticity(vInfo, *grid, uBody, oldAccVort, rhoS);
			
			if (step>50)
			{
				timer.start();
				vOld = uBody[1];
				computeBodyVelocity(vInfo, *grid, uBody, omegaBody, rhoS, gravity, dt, lambda);
				shape->updatePosition(uBody, omegaBody, dt);
				processOMP<OperatorComputeShape>(shape, vInfo, *grid);
				timeUB += timer.stop();
			}
			
			// penalization
			timer.start();
			Real centerOfMass[2] = {0,0};
			shape->getPosition(centerOfMass);
			processOMP<OperatorPenalization>(dt,uBody[0],uBody[1],omegaBody,centerOfMass[0],centerOfMass[1],lambda,vInfo,*grid);
			timePenalization += timer.stop();
			
			// body integrals
			timer.start();
			timeUB += timer.stop();
			
			if (step<100)
			{
				// this still needs to be corrected to the frame of reference!
				double accM = (uBody[1]-vOld)/dt;
				double accT = (rhoS-1)/(rhoS+1) * gravity[1];
				double accN = (rhoS-1)/(rhoS) * gravity[1];
				cout << "Acceleration with added mass (measured, expected, no added mass)\t" << accM << "\t" << accT << "\t" << accN << endl;
				stringstream ss;
				ss << path2file << "_addedmass.dat";
				ofstream myfile(ss.str(), fstream::app);
				myfile << step << " " << time << " " << accM << " " << accT << " " << accN << endl;
			}
			
			// compute diagnostics
			if (step % 10 == 0)
			{
				timer.start();
				_diagnostics();
				timeDiagnostics += timer.stop();
			}
			
			//dump some time steps every now and then
			if (time>nextDumpTime)
			//if (step % 1 == 0)
			{
				nextDumpTime += .05;
				
				timer.start();
				stringstream ss;
				ss << path2file << "-" << step << ".vti";
				cout << ss.str() << endl;
				
				dumper.Write(*grid, ss.str());
				
				Layer vorticity(sizeX,sizeY,1);
				processOMP<Lab, OperatorVorticity>(vorticity,vInfo,*grid);
				stringstream sVort;
				sVort << path2file << "Vorticity-" << step << ".vti";
				dumpLayer2VTK(step,sVort.str(),vorticity,1);
				timeDump += timer.stop();
			}
			
			if (step % 100 == 0)
			{
				double totalTime = timeDT + timeAdvection + timeDiffusion + timePressure + timePenalization + timeDiagnostics + timeDump;
				cout << "=== Timing Report ===\n";
				cout << "\tDT\t\t" << setprecision(3) << timeDT << "s ( " << setprecision(2) << 100.*timeDT/totalTime << " % )\n";
				cout << "\tAdvection\t" << setprecision(3) << timeAdvection << "s ( " << setprecision(2) << 100.*timeAdvection/totalTime << " % )\n";
				cout << "\tDiffusion\t" << setprecision(3) << timeDiffusion << "s ( " << setprecision(2) << 100.*timeDiffusion/totalTime << " % )\n";
				cout << "\tPressure\t" << setprecision(3) << timePressure << "s ( " << setprecision(2) << 100.*timePressure/totalTime << " % )\n";
				//cout << "\t\tDivergence\t" << timeDivergence << endl;
				//cout << "\t\tCut\t\t" << timeCut << endl;
				//cout << "\t\tPoisson\t\t" << timePoisson << endl;
				//cout << "\t\tUpdate\t\t" << timeUpdatePressure << endl;
				cout << "\tPenalization\t" << setprecision(3) << timePenalization << "s ( " << setprecision(2) << 100.*timePenalization/totalTime << " % )\n";
				cout << "\tDiagnostics\t" << setprecision(3) << timeDiagnostics << "s ( " << setprecision(2) << 100.*timeDiagnostics/totalTime << " % )\n";
				cout << "\tDump\t\t" << setprecision(3) << timeDump << "s ( " << setprecision(2) << 100.*timeDump/totalTime << " % )\n";
			}
			
			time += dt;
			
			// check nondimensional time
			if (time>tEnd)
			{
				timer.start();
				stringstream ss;
				ss << path2file << "-Final.vti";
				cout << ss.str() << endl;
				
				dumper.Write(*grid, ss.str());
				
				Layer vorticity(sizeX,sizeY,1);
				processOMP<Lab, OperatorVorticity>(vorticity,vInfo,*grid);
				stringstream sVort;
				sVort << path2file << "Vorticity-Final.vti";
				dumpLayer2VTK(step,sVort.str(),vorticity,1);
				timeDump += timer.stop();
				
				
				
				double totalTime = timeDT + timeAdvection + timeDiffusion + timePressure + timePenalization + timeDiagnostics + timeDump;
				cout << "=== Final Timing Report ===\n";
				cout << "\tDT\t\t" << setprecision(3) << timeDT << "s ( " << setprecision(2) << 100.*timeDT/totalTime << " % )\n";
				cout << "\tAdvection\t" << setprecision(3) << timeAdvection << "s ( " << setprecision(2) << 100.*timeAdvection/totalTime << " % )\n";
				cout << "\tDiffusion\t" << setprecision(3) << timeDiffusion << "s ( " << setprecision(2) << 100.*timeDiffusion/totalTime << " % )\n";
				cout << "\tPressure\t" << setprecision(3) << timePressure << "s ( " << setprecision(2) << 100.*timePressure/totalTime << " % )\n";
				//cout << "\t\tDivergence\t" << timeDivergence << endl;
				//cout << "\t\tCut\t\t" << timeCut << endl;
				//cout << "\t\tPoisson\t\t" << timePoisson << endl;
				//cout << "\t\tUpdate\t\t" << timeUpdatePressure << endl;
				cout << "\tPenalization\t" << setprecision(3) << timePenalization << "s ( " << setprecision(2) << 100.*timePenalization/totalTime << " % )\n";
				cout << "\tDiagnostics\t" << setprecision(3) << timeDiagnostics << "s ( " << setprecision(2) << 100.*timeDiagnostics/totalTime << " % )\n";
				cout << "\tDump\t\t" << setprecision(3) << timeDump << "s ( " << setprecision(2) << 100.*timeDump/totalTime << " % )\n";
				
				break;
			}
		}
	}
}