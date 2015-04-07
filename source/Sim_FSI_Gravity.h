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

class Sim_FSI_Gravity : public Simulation_FSI
{
protected:
	Real uBody[2], omegaBody;
	double dtCFL, dtFourier, dtBody;
	double re, nu;
	double minRho;
	bool bSplit;
	
	Real gravity[2];
	
	int stepStartBody; // could become time
	
	void _diagnostics();
	void _ic();
	double _nonDimensionalTime();
	
	void _outputSettings(ostream& outStream);
	void _inputSettings(istream& inStream);
	
	// should this stuff be moved? - serialize method will do that
	void _dumpSettings(ostream& outStream);
	
public:
	Sim_FSI_Gravity(const int argc, const char ** argv);
	~Sim_FSI_Gravity();
	
	void init();
	void simulate();
};

#endif /* defined(__CubismUP_2D__Sim_FSI_Gravity__) */
