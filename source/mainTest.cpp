//
//  mainTest.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
using namespace std;

#include "common.h"

#include "TestDiffusion.h"
#include "TestAdvection.h"
#include "TestPressure.h"
#include "TestVarCoeffPoisson.h"
#include "TestGravity.h"
#include "TestPenalization.h"
#include "TestTranslation.h"
#include "TestRotation.h"
#include "TestTravelingWave.h"
#include "TestShearLayer.h"
#include "Definitions.h"

int main(int argc, const char **argv)
{
#ifdef _MULTIGRID_
	MPI_Init(&argc, &argv);
#endif // _MULTIGRID_
	
	ArgumentParser parser(argc,argv);
	int solver = parser("-solver").asInt(0);
	int ic = parser("-ic").asInt(0);
	parser.set_strict_mode();
	
	string test = parser("-test").asString();
	const int minBPD = parser("-minBPD").asInt();
	const int maxBPD = parser("-maxBPD").asInt();
	
	if (test=="advection")
	{
		if (ic==0)
		{
			cout << "========================================================================================\n";
			cout << "\t\tAdvection Test - Linear Field\n";
			cout << "========================================================================================\n";
		}
		else if (ic==1)
		{
			cout << "========================================================================================\n";
			cout << "\t\tAdvection Test - Vortex Field\n";
			cout << "========================================================================================\n";
		}
		else
		{
			cout << "IC " << ic << " doesn't exist\n";
			abort();
		}
		
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestAdvection * advection = new TestAdvection(argc, argv, ic, bpd);
			advection->run();
			advection->check();
			delete advection;
		}
	}
	else if (test=="diffusion")
	{
		cout << "========================================================================================\n";
		cout << "\t\tDiffusion Test\n";
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestDiffusion * diffusion = new TestDiffusion(argc, argv, bpd);
			diffusion->run();
			diffusion->check();
			delete diffusion;
		}
	}
	else if (test=="pressure")
	{
		cout << "========================================================================================\n";
		if (solver==0)
			cout << "\t\tPressure Test - Stencil - ";
		else if (solver==1)
			cout << "\t\tPressure Test - Spectral - ";
#ifdef _MULTIGRID_
		else if (solver==2)
			cout << "\t\tPressure Test - Multigrid (Constant Coefficents) - ";
		else if (solver==3)
			cout << "\t\tPressure Test - Multigrid (Variable Coefficents) - ";
#endif // _MULTIGRID_
		else
		{
			cout << "Solver case " << solver << " doesn't exist\n";
			abort();
		}
		
		if (ic==0)
			cout << "Poisson equation\n";
		else if (ic==1)
			cout << "Velocity field (Projection)\n";
		else if (ic==2)
			cout << "Mixed Boundary Conditions\n";
		else
			abort();
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestPressure * pressure = new TestPressure(argc, argv, solver, ic, bpd);
			pressure->run();
			pressure->check();
			delete pressure;
		}
	}
	else if (test=="poisson")
	{
		cout << "========================================================================================\n";
		cout << "\t\tPoisson Test - Multigrid - ";
		
		if (ic==0)
			cout << "Constant Coefficients\n";
		else if (ic==1)
			cout << "Variable Coefficients\n";
		else
			abort();
		cout << "========================================================================================\n";
		
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestVarCoeffPoisson * poisson = new TestVarCoeffPoisson(argc, argv, ic, bpd);
			poisson->run();
			poisson->check();
			delete poisson;
		}
	}
	else if (test=="gravity")
	{
		cout << "========================================================================================\n";
		cout << "\t\tGravity Test\n";
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestGravity * gravity = new TestGravity(argc, argv, bpd);
			gravity->run();
			gravity->check();
			delete gravity;
		}
	}
	else if (test=="penalization")
	{
		cout << "========================================================================================\n";
		cout << "\t\tPenalization Test\n";
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestPenalization * penalization = new TestPenalization(argc, argv, bpd);
			penalization->run();
			delete penalization;
		}
	}
	else if (test=="translation")
	{
		cout << "========================================================================================\n";
		if (ic==0)
		cout << "\t\tTranslation Test - Constant Velocity\n";
		else if (ic==1)
		cout << "\t\tTranslation Test - Velocity from Flow\n";
		else
		abort();
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestTranslation * translation = new TestTranslation(argc, argv, ic, bpd);
			translation->run();
			translation->check();
			delete translation;
		}
	}
	else if (test=="rotation")
	{
		cout << "========================================================================================\n";
		if (ic==0)
		cout << "\t\tRotation Test - Constant Velocity\n";
		else if (ic==1)
		cout << "\t\tRotation Test - Velocity from Flow\n";
		else
		abort();
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestRotation * rotation = new TestRotation(argc, argv, ic, bpd);
			rotation->run();
			rotation->check();
			delete rotation;
		}
	}
	else if (test=="travelingwave")
	{
		cout << "========================================================================================\n";
		cout << "\t\tTraveling Wave Test\n";
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestTravelingWave * wave = new TestTravelingWave(argc, argv, bpd);
			wave->run();
			wave->check();
			delete wave;
		}
	}
	else if (test=="shearlayer")
	{
		cout << "========================================================================================\n";
		cout << "\t\tShear Layer Test\n";
		cout << "========================================================================================\n";
		/*
		// reference run
		TestShearLayer * shearlayer = new TestShearLayer(argc, argv, 32);
		shearlayer->run();
		delete shearlayer;
		//*/
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestShearLayer * shearlayer = new TestShearLayer(argc, argv, bpd);
			shearlayer->run();
			shearlayer->check();
			delete shearlayer;
		}
	}
	
	
#ifdef _MULTIGRID_
	MPI_Finalize();
#endif // _MULTIGRID_
	
    return 0;
}