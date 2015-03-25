//
//  Sim_FSI_Fixed.h
//  CubismUP_2D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_2D__Sim_FSI_Fixed__
#define __CubismUP_2D__Sim_FSI_Fixed__

#include "Simulation_FSI.h"

class Sim_FSI_Fixed : public Simulation_FSI
{
protected:
    double uinf, re, nu;
	double dtCFL, dtFourier;
	
	void _diagnostics();
	void _ic();
	double _nonDimensionalTime();
	
public:
	Sim_FSI_Fixed(const int argc, const char ** argv);
	
	virtual ~Sim_FSI_Fixed();
	
    void simulate();
};


#endif /* defined(__CubismUP_2D__Sim_FSI_Fixed__) */
