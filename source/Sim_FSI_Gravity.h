//
//  Sim_FSI_Gravity.h
//  CubismUP_2D
//
//  Created by Christian Conti on 1/26/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_2D__Sim_FSI_Gravity__
#define __CubismUP_2D__Sim_FSI_Gravity__

#include "Simulation_FSI.h"
#ifdef _MULTIGRID_
#include "MultigridHypre.h"
#endif

class Sim_FSI_Gravity : public Simulation_FSI
{
protected:
	Real uBody[2], omegaBody;
	double dtCFL, dtFourier, dtBody;
	double re, nu;
	double CFL;
	double minRho;
	bool bSplit;
	
	Real gravity[2];
	
#ifdef _MULTIGRID_
	MultigridHypre mg;
	int rank, nprocs;
#endif
	
	void _diagnostics();
	void _ic();
	double _nonDimensionalTime();
	
	// should this stuff be moved?
	void _dumpSettings(ostream& mystream);
	void _dumpDivergence(const int step, const Real rho0, const Real dt);
	void _solvePressure();
	
public:
	Sim_FSI_Gravity(const int argc, const char ** argv);
	
	~Sim_FSI_Gravity();

	void simulate();
};

#endif /* defined(__CubismUP_2D__Sim_FSI_Gravity__) */
