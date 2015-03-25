//
//  Sim_FSI_Fixed.h
//  CubismUP_2D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_2D__Sim_FSI_Fixed__
#define __CubismUP_2D__Sim_FSI_Fixed__

#include <stdio.h>
#include "Definitions.h"
#include "Timer.h"
#include "ProcessOperatorsOMP.h"
#include "OperatorIC.h"
#include "OperatorAdvection.h"
#include "OperatorDiffusion.h"
#include "OperatorPenalization.h"
#include "OperatorDivergence.h"
#include "OperatorVorticity.h"
#include "Operators_DFT.h"
#include "PoissonSolverScalarFFTW.h"
#include "OperatorGradP.h"
#include "OperatorComputeShape.h"

class Sim_FSI_Fixed
{
protected:
    ArgumentParser parser;
    
    double uinf, re, nu;
	double lambda;
	double dt, dtCFL, dtFourier, time, tEnd;
	double rhoS;
	Shape * shape;
	
    int step, nsteps;
    int bpdx, bpdy;
    
    string path2file;
    SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;

    FluidGrid * grid;
    
    // add states
	// simulate should only be called after setup
	
	virtual void _diagnostics();
	double _nonDimensionalTime(double u);
	
public:
    Sim_FSI_Fixed(const int argc, const char ** argv) : parser(argc, argv), uinf(0), step(-1), time(-1), nsteps(0), re(0), nu(0), lambda(0)
	{
		int rank, nprocs;
#ifdef _MULTIGRID_
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
#endif
		
		if (rank==0)
		{
			cout << "====================================================================================================================\n";
			cout << "\t\t\tFlow past a fixed cylinder\n";
			cout << "====================================================================================================================\n";
		}
    }
	
    virtual ~Sim_FSI_Fixed()
	{
		delete grid;
		delete shape;
	}
    
    virtual void setup();
    virtual void simulate();
};


#endif /* defined(__CubismUP_2D__Sim_FSI_Fixed__) */
