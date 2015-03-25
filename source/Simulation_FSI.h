//
//  Simulation_FSI.h
//  CubismUP_2D
//
//	Base class for Fluid-Structure Interaction (FSI) simulations from which any FSI simulation case should inherit
//	Contains the base structure and interface that any FSI simulation class should have
//	Inherits from Simulation_Fluid
//	Assumes use of Penalization to handle rigid body
//
//  Created by Christian Conti on 3/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_Simulation_FSI_h
#define CubismUP_2D_Simulation_FSI_h

#include "Simulation_Fluid.h"
#include "Shape.h"
#include "Operators_DFT.h"

class Simulation_FSI : public Simulation_Fluid
{
protected:
	// penalization parameter
	double lambda;
	
	// body
	Shape * shape;
	
public:
	Simulation_FSI(const int argc, const char ** argv) : Simulation_Fluid(argc,argv)
	{
		lambda = parser("-lambda").asDouble(1e8);
		
		double rhoS = parser("-rhoS").asDouble(1);
		Real centerOfMass[2] = {0,0};
		bool bPeriodic[2] = {false,false};
		
		string shapeType = parser("-shape").asString("disk");
		if (shapeType=="disk")
		{
			Real radius = parser("-radius").asDouble(0.1);
			shape = new Disk(centerOfMass, radius, rhoS, 2, 2, bPeriodic);
		}
		else if (shapeType=="ellipse")
		{
			Real semiAxis[2] = {parser("-semiAxisX").asDouble(0.1),parser("-semiAxisY").asDouble(0.2)};
			Real angle = parser("-angle").asDouble(0.0);
			shape = new Ellipse(centerOfMass, semiAxis, angle, rhoS, 2, 2, bPeriodic);
		}
		else
		{
			cout << "Error - this shape is not currently implemented! Aborting now\n";
			abort();
		}
	}
	
	virtual ~Simulation_FSI()
	{
		delete shape;
	}
};

#endif
