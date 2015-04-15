//
//  Sim_RayleighTaylor.h
//  CubismUP_2D
//
//  Created by Christian Conti on 4/10/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_2D__Sim_RayleighTaylor__
#define __CubismUP_2D__Sim_RayleighTaylor__

#include "Simulation_MP.h"

class Sim_RayleighTaylor : public Simulation_MP
{
protected:
	void _diagnostics();
	void _ic();
	double _nonDimensionalTime();
	
	void _outputSettings(ostream& outStream);
	void _inputSettings(istream& inStream);
	
public:
	Sim_RayleighTaylor(const int argc, const char ** argv);
	~Sim_RayleighTaylor();
	
	void init();
	void simulate();
};

#endif /* defined(__CubismUP_2D__Sim_RayleighTaylor__) */
