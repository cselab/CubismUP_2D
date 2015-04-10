//
//  Sim_Multiphase.h
//  CubismUP_2D
//
//  Created by Christian Conti on 4/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_2D__Sim_Multiphase__
#define __CubismUP_2D__Sim_Multiphase__

#include "Simulation_Fluid.h"

class Sim_Multiphase : public Simulation_Fluid
{
protected:
	double dtCFL, dtFourier;
	double re, nu;
	double minRho, rhoS;
	bool bSplit;
	
	Real gravity[2];
	
	void _diagnostics();
	void _ic();
	double _nonDimensionalTime();
	
	void _outputSettings(ostream& outStream);
	void _inputSettings(istream& inStream);
	
	// should this stuff be moved? - serialize method will do that
	void _dumpSettings(ostream& outStream);
	
public:
	Sim_Multiphase(const int argc, const char ** argv);
	~Sim_Multiphase();
	
	void init();
	void simulate();
};

#endif /* defined(__CubismUP_2D__Sim_Multiphase__) */
