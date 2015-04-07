//
//  main.cpp
//  CubismUP_2D
//
//  Created by Christian Conti on 1/7/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
using namespace std;

#include "Definitions.h"

#include "Simulation_Fluid.h"
#include "Sim_FSI_Fixed.h"
#include "Sim_FSI_Moving.h"
#include "Sim_FSI_Gravity.h"

int main(int argc, const char **argv)
{
#ifdef _MULTIGRID_
	MPI_Init(&argc, &argv);
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	if (rank==0)
#endif // _MULTIGRID_
	{
		cout << "====================================================================================================================\n";
		cout << "\t\tCubism UP 2D (velocity-pressure 2D incompressible Navier-Stokes solver)\n";
		cout << "====================================================================================================================\n";
	}
	
	ArgumentParser parser(argc,argv);
	parser.set_strict_mode();
	
	string simSetting = parser("-sim").asString();
	Simulation_Fluid * sim;
	
	if (simSetting=="fixed")
		sim = new Sim_FSI_Fixed(argc, argv);
	else if (simSetting=="moving")
		sim = new Sim_FSI_Moving(argc, argv);
	else if (simSetting=="falling")
		sim = new Sim_FSI_Gravity(argc, argv);
	
	sim->init();
	sim->simulate();
	
#ifdef _MULTIGRID_
	MPI_Finalize();
#endif // _MULTIGRID_
	
	return 0;
}