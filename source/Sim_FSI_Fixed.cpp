//
//  Sim_FSI_Fixed.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Sim_FSI_Fixed.h"


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
	//cout << step << " " << time << " " << bpd << " " << dt << " " << dtCFL << " " << dtFourier << " " << drag << " " << lambda << endl;
	myfile << step << " " << _nonDimensionalTime(uinf) << " " << bpdx << " " << dt << " " << dtCFL << " " << dtFourier << " " << cD << " " << lambda << endl;
}

double Sim_FSI_Fixed::_nonDimensionalTime(double u)
{
	return time*abs(u)/shape->getCharLength();
}

void Sim_FSI_Fixed::setup()
{
	// simulation settings
	nsteps = parser("-nsteps").asInt(1e9);
	uinf = parser("-uinf").asDouble(0.1);
	re = parser("-Re").asDouble(100);
	lambda = parser("-lambda").asDouble(1e4);
	tEnd = parser("-tend").asDouble(50);
	rhoS = parser("-rhoS").asDouble(1);
	
	// output settings
    path2file = parser("-file").asString("../data/FlowPastCylinder");
    
    // computational space settings
	parser.set_strict_mode();
	bpdx = parser("-bpdx").asInt();
	bpdy = parser("-bpdy").asInt();
	
	if (bpdx!=bpdy)
	{
		cout << "Assuming square domain for periodicity in OperatorComputeDisk\n";
		abort();
	}
	
	Real center[2] = {.15,.5};
	Real radius = parser("-radius").asDouble(0.1);
	bool bPeriodic[2] = {false,false};
	shape = new Disk(center, radius, rhoS, 2, 2, bPeriodic);
	nu = shape->getCharLength()*uinf/re;
	
    grid = new FluidGrid(bpdx,bpdy,1);
	assert(grid != NULL);
	
	Timer timerIC;
	
	// setup initial conditions
	timerIC.start();
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
	OperatorIC ic(shape, uinf);
	FluidBlockProcessing::process(vInfo, ic, grid);
	double timeIC = timerIC.stop();
	
    stringstream ss;
    ss << path2file << "-IC.vti" ;
    cout << ss.str() << endl;
    
    dumper.Write(*grid, ss.str());
	cout << "Time IC (u, p):\t" << timeIC << endl;// << "\t" << timePoissonIC << endl;
}

void Sim_FSI_Fixed::simulate()
{
    time = 0;
	Real uBody[2] = {0,0};
	double maxU = uinf;
	Timer timer, timerPressure;
	
	double timeDT = 0;
	double timeAdvection = 0;
	double timeDiffusion = 0;
	double timePressure = 0;
	double timePenalization = 0;
	double timeDiagnostics = 0;
	double timeDump = 0;
	
	double timeDivergence = 0;
	double timePoisson = 0;
	double timeUpdatePressure = 0;
	
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
#ifdef _SP_COMP_
	PoissonSolverScalarFFTW<FluidGrid, StreamerDiv> pressureSolver(NTHREADS);
#else
	cout << "FFTW double precision not supported - aborting now!\n";
	abort();
#endif
	
    for(step=0; step<nsteps; ++step)
	{
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		
        // choose dt (CFL, Fourier)
		timer.start();
		maxU = findMaxUOMP(vInfo,*grid);
		dtFourier = .1*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
		dtCFL     = .1*vInfo[0].h_gridpoint/abs(maxU);
        dt = min(dtCFL,dtFourier);
		//cout << "dt (Fourier, CFL): " << dtFourier << " " << dtCFL << endl;
		timeDT += timer.stop();
        
        // the operators could be put into a vector pipeline
		
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
		
		// pressure
		timer.start();
		timerPressure.start();
		processOMP<Lab, OperatorDivergence>(dt, vInfo, *grid);
		timeDivergence += timerPressure.stop();
		
		timerPressure.start();
#ifdef _SP_COMP_
		pressureSolver.solve(*grid,true);
#else
		cout << "FFTW double precision not supported - aborting now!\n";
		abort();
#endif
		timePoisson += timerPressure.stop();
		
		timer.start();
		processOMP<Lab, OperatorGradP>(dt, vInfo, *grid);
		timeUpdatePressure += timerPressure.stop();
		timePressure += timer.stop();
		 
		// penalization
		timer.start();
		Real g[2] = {0,0};
		Real centerOfMass[2] = {0,0};
		shape->getPosition(centerOfMass);
		processOMP<OperatorPenalization>(dt,0,0,0,centerOfMass[0],centerOfMass[1],lambda,vInfo,*grid);
		timePenalization += timer.stop();
		
		// compute diagnostics
		if (step % 10 == 0)
		{
			timer.start();
			_diagnostics();
			timeDiagnostics += timer.stop();
		}
        
        //dump some time steps every now and then
        if(step % 1000 == 0)
        {
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
		if (_nonDimensionalTime(uinf)>tEnd)
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