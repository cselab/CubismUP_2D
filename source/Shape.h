//
//  Shape.h
//  CubismUP_2D
//
//	Virtual shape class which defines the interface
//	Default simple geometries are also provided and can be used as references
//
//	This class only contains static information (position, orientation,...), no dynamics are included (e.g. velocities,...)
//
//  Created by Christian Conti on 3/6/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_Shape_h
#define CubismUP_2D_Shape_h

class Shape
{
protected:
	Real center[2], orientation;
	const Real rhoS;
	
	const Real mollChi;
	const Real mollRho; // currently not used - need to change in rho method
	
	Real smoothHeaviside(Real rR, Real radius, Real eps)
	{
		if (rR < radius-eps*.5)
			return (Real) 1.;
		else if (rR > radius+eps*.5)
			return (Real) 0.;
		else
			return (Real) ((1.+cos(M_PI*((rR-radius)/eps+.5)))*.5);
	}
	
	// TODO: need to include boundary conditions in case of periodicity
	
public:
	Shape(Real center[2], Real orientation, const Real rhoS, const Real mollChi, const Real mollRho, bool bPeriodic[2]) : center{center[0],center[1]}, orientation(orientation), rhoS(rhoS), mollChi(mollChi), mollRho(mollRho)
	{
		if (bPeriodic[0] || bPeriodic[1])
		{
			cout << "Periodic shapes are currently unsupported\n";
			abort();
		}
	}
	
	virtual ~Shape() {}
	
	virtual Real chi(Real p[2], Real h) const = 0;
	virtual Real getCharLength() const = 0;
	
	
	void updatePosition(Real u[2], Real omega, Real dt)
	{
		center[0] += dt*u[0];
#ifndef _MOVING_FRAME_
        center[1] += dt*u[1];
#endif
		orientation += dt*omega;
	}
	
	void setPosition(Real com[2])
	{
		center[0] = com[0];
		center[1] = com[1];
	}
	
	void getPosition(Real com[2])
	{
		com[0] = center[0];
		com[1] = center[1];
	}
	
	Real getOrientation() const
	{
		return orientation;
	}
	
	inline Real getRhoS() const
	{
		return rhoS;
	}
	
	Real rho(Real p[2], Real h)
	{
		Real mask = chi(p,h);
		
		return rhoS*mask + 1.*(1.-mask);
	}
	
	virtual void outputSettings(ostream &outStream)
	{
		outStream << "centerX " << center[0] << endl;
		outStream << "centerY " << center[1] << endl;
		outStream << "orientation " << orientation << endl;
		outStream << "rhoS " << rhoS << endl;
		outStream << "mollChi " << mollChi << endl;
		outStream << "mollRho " << mollRho << endl;
	}
};

class Disk : public Shape
{
protected:
	Real radius;
	
public:
	Disk(Real center[2], Real radius, const Real rhoS, const Real mollChi, const Real mollRho, bool bPeriodic[2]) : Shape(center, 0, rhoS, mollChi, mollRho, bPeriodic), radius(radius) {}
	
	Real chi(Real p[2], Real h) const
	{
		const Real d[2] = { abs(p[0]-center[0]), abs(p[1]-center[1]) };
		const Real dist = sqrt(d[0]*d[0]+d[1]*d[1]);
		
		return smoothHeaviside(dist, radius, mollChi*sqrt(2)*h);
	}
	
	Real getCharLength() const
	{
		return 2 * radius;
	}
	
	void outputSettings(ostream &outStream)
	{
		outStream << "Disk\n";
		outStream << "radius " << radius << endl;
		
		Shape::outputSettings(outStream);
	}
};

class Ellipse : public Shape
{
protected:
	// these quantities are defined in the local coordinates of the ellipse
	Real semiAxis[2];
	
public:
	Ellipse(Real center[2], Real semiAxis[2], Real orientation, const Real rhoS, const Real mollChi, const Real mollRho, bool bPeriodic[2]) : Shape(center, orientation, rhoS, mollChi, mollRho, bPeriodic), semiAxis{semiAxis[0],semiAxis[1]} {}
	
	Real chi(Real p[2], Real h) const
	{
		const Real angle = atan2(p[1]-center[1],p[0]-center[0]) - orientation;
		const Real x = semiAxis[0]*cos(angle);
		const Real y = semiAxis[1]*sin(angle);
		const Real radius = semiAxis[0]*semiAxis[1] / sqrt(x*x + y*y);
		const Real dist = sqrt( (p[0]-center[0])*(p[0]-center[0]) + (p[1]-center[1])*(p[1]-center[1]) );
		return smoothHeaviside(dist, radius, mollChi*sqrt(2)*h);
	}
	
	Real getCharLength() const
	{
		return 2 * semiAxis[1];
	}
	
	void outputSettings(ostream &outStream)
	{
		outStream << "Ellipse\n";
		outStream << "semiAxisX " << semiAxis[0] << endl;
		outStream << "semiAxisY " << semiAxis[1] << endl;
		
		Shape::outputSettings(outStream);
	}
};

#endif
