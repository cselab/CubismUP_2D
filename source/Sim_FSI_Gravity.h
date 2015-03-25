//
//  Sim_FSI_Gravity.h
//  CubismUP_2D
//
//  Created by Christian Conti on 1/26/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_2D__Sim_FSI_Gravity__
#define __CubismUP_2D__Sim_FSI_Gravity__

#include <stdio.h>
#include "Definitions.h"
#include "Sim_FSI_Fixed.h"
#include "OperatorGravity.h"
#include "OperatorSplitP.h"
#ifdef _MULTIGRID_
#include "MultigridHypre.h"
#endif

class Sim_FSI_Gravity : public Sim_FSI_Fixed
{
protected:
	Real uBody[2];
	Real omegaBody;
	Real gravity[2];
	
	double dtBody;
	double minRho;
	double CFL;
	
#ifdef _MULTIGRID_
	MultigridHypre mg;
	int rank, nprocs;
	bool bSplit;
#endif
	
	void _diagnostics();
	void _dumpSettings(ostream& mystream);
	void _dumpDivergence(const int step, const Real rho0, const Real dt);
	void _solvePressure();
	
public:
	Sim_FSI_Gravity(const int argc, const char ** argv) : Sim_FSI_Fixed(argc, argv), uBody{0,0}, omegaBody(0), gravity{0,-9.81}
	{
#ifdef _MULTIGRID_
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
		bSplit = parser("-split").asBool(false);
		
		if (rank==0)
#endif
		{
			cout << "====================================================================================================================\n";
			cout << "\t\t\tFlow past a falling cylinder\n";
			cout << "====================================================================================================================\n";
		}
	}
	
	~Sim_FSI_Gravity() {}
	
	void setup();
	void simulate();
};

#endif /* defined(__CubismUP_2D__Sim_FSI_Gravity__) */
