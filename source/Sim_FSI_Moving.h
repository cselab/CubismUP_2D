//
//  Sim_FSI_Moving.h
//  CubismUP_2D
//
//  Created by Christian Conti on 1/26/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_2D__Sim_FSI_Moving__
#define __CubismUP_2D__Sim_FSI_Moving__

#include <stdio.h>
#include "Definitions.h"
#include "Sim_FSI_Fixed.h"

class Sim_FSI_Moving : public Sim_FSI_Fixed
{
protected:
	Real uBody[2];
	double dtBody;
	
	void _diagnostics();
	
public:
	Sim_FSI_Moving(const int argc, const char ** argv) : Sim_FSI_Fixed(argc, argv), uBody{0,0}
	{
		cout << "====================================================================================================================\n";
		cout << "\t\t\tFlow past a moving cylinder\n";
		cout << "====================================================================================================================\n";
	}
	
	~Sim_FSI_Moving() {}
	
	void setup();
	void simulate();
};

#endif /* defined(__CubismUP_2D__Sim_FSI_Moving__) */
