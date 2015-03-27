/*
 *  Operators_DFT.h
 *  VM2D
 *
 *  Created by Diego Rossinelli on 2/10/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

#define _USE_MATH_DEFINES

#include <math.h>
#include <typeinfo>
#include <xmmintrin.h>
#include <emmintrin.h>

#include "Layer.h"
#include "Matrix2D.h"

template<typename T> struct FFTW_TYPES {};

#include <fftw3.h>

template<> struct FFTW_TYPES<float> { 
	typedef fftwf_complex COMPLEX; 
	typedef fftwf_plan PLAN;
	
	struct ALLOCATOR
	{
		template<typename T>
		static void alloc(T*& ptr, const int N) { ptr = (T *)fftwf_malloc(N*sizeof(T));}
		template<typename T>
		static void free(T*& ptr) { fftwf_free(ptr);}
	};
	
	struct PLANNER
	{
		static void del(PLAN plan) { fftwf_destroy_plan(plan); }
		static void ex(PLAN plan) { fftwf_execute(plan); }
		static PLAN c2r(const int nX, const int nY, COMPLEX* in, float * out) {return fftwf_plan_dft_c2r_2d(nY, nX, in, out, FFTW_MEASURE);}
		static PLAN r2c(const int nX, const int nY, float * in, COMPLEX * out){ return fftwf_plan_dft_r2c_2d(nY, nX, in, out, FFTW_MEASURE);}
	};
};

template<> struct FFTW_TYPES<double> {
	typedef fftw_complex COMPLEX; 
	typedef fftw_plan PLAN;
	
	struct ALLOCATOR
	{
		template<typename T>
		static void alloc(T*& ptr, const int N) {  ptr = (T*)_mm_malloc(sizeof(T)*N, 16); assert(((unsigned long int)ptr & 0xf) == 0);/* sucks ptr = (T *)fftw_malloc(N*sizeof(T)); */}
		template<typename T>
		static void free(T*& ptr) { _mm_free(ptr); /*fftw_free(ptr);*/}
	};
	
	struct PLANNER
	{
		static void del(PLAN plan) { fftw_destroy_plan(plan); }
		static void ex(PLAN plan) { fftw_execute(plan); }
		static PLAN c2r(const int nX, const int nY, COMPLEX* in, double * out) {return fftw_plan_dft_c2r_2d(nY, nX, in, out, FFTW_MEASURE);}
		static PLAN r2c(const int nX, const int nY, double * in, COMPLEX * out){ return fftw_plan_dft_r2c_2d(nY, nX, in, out, FFTW_MEASURE);}
	};	
};

template<typename T> struct FFTW_INIT {};
template<> struct FFTW_INIT<float>  {	FFTW_INIT() {fftwf_init_threads(); fftwf_plan_with_nthreads(NTHREADS); }; };
template<> struct FFTW_INIT<double>  {	FFTW_INIT() {fftw_init_threads(); fftw_plan_with_nthreads(NTHREADS); printf("init fftw with %d threads\n", NTHREADS); }; };


class VelocitySolver_FFTW
{
private:
	//forbidden
	VelocitySolver_FFTW(const VelocitySolver_FFTW&):psi_physical(NULL), fftw_metadata(NULL), m_bDoublePrecision(false)
	{
		abort();
	}
	
	VelocitySolver_FFTW& operator=(VelocitySolver_FFTW& ){abort(); return *this;}
	
protected:
	bool m_bDoublePrecision;
	int m_nSize[2];
	Real (* data_fourier)[2];
	Real * psi_physical;
	void * fftw_metadata;
	
	template<typename RealType> 
	void _init(RealType * vorticity, RealType (*&data_fourier)[2], RealType * psi_physical );
	
	void _solve_phi(int sizeX, int sizeY, Real H);
	
	void _solve_phi_core(Real (*&data_fourier)[2], Real H);
	
	void _dispose();
	
	template <int FDdir>
	void _solve_velocity(const Layer& psi, Layer& velocity)
	{
		const double invH2[2] = {
			1./psi.getH0(),
			1./psi.getH1()
		};
		
		const Real factor[2] = {
			invH2[0]*0.5,
			invH2[1]*0.5
		};
		
		if (FDdir==0)
		{
			for(int dy=0; dy<psi.sizeY; dy++)
			for(int dx=0; dx<psi.sizeX; dx++)
			{
				const Real psi_right = psi.read((dx + 1) % psi.sizeX, dy);
				const Real psi_left = psi.read((dx - 1 + psi.sizeX) % psi.sizeX, dy);
				const Real psi_top = psi.read(dx, (dy + 1) % psi.sizeY);
				const Real psi_bottom = psi.read(dx, (dy - 1 + psi.sizeY) % psi.sizeY);
				//const Real psi_center =  psi.read(dx, dy);
				
				velocity(dx, dy, 0) =  factor[1]*(psi_top - psi_bottom);
				velocity(dx, dy, 1) =  -factor[0]*(psi_right - psi_left);
			}
		}
		else if (FDdir==1)
		{
			for(int dy=0; dy<psi.sizeY; dy++)
				for(int dx=0; dx<psi.sizeX; dx++)
				{
					const Real psi_right = psi.read((dx + 1) % psi.sizeX, dy);
					const Real psi_top = psi.read(dx, (dy + 1) % psi.sizeY);
					const Real psi_center =  psi.read(dx, dy);
					
					velocity(dx, dy, 0) =  invH2[1]*(psi_top - psi_center);
					velocity(dx, dy, 1) =  -invH2[0]*(psi_right - psi_center);
				}
		}
		else
			abort();
			
	}
	
public:
	VelocitySolver_FFTW(bool bDoublePrecision, Layer& vorticity);
	
	~VelocitySolver_FFTW()
	{
		_dispose();
	}
	
	void operator() (Layer& velocity, int FiniteDifferenceDir=0)
	{
		assert(velocity.sizeX == m_nSize[0]);
		assert(velocity.sizeY == m_nSize[1]);
		
		assert(velocity.getH0() == velocity.getH1());
		_solve_phi(velocity.sizeX, velocity.sizeY, velocity.getH0());
		
		const Layer& psi_as_layer = *((const Layer*)psi_physical);
		
		if (FiniteDifferenceDir==0)
			_solve_velocity<0>(psi_as_layer, velocity);
		else if (FiniteDifferenceDir == 1)
			_solve_velocity<1>(psi_as_layer, velocity);
		else
			abort();
		
	}
};

class PressureSolver_FFTW
{
protected:
	int m_nSize[2];
	bool m_bReady;
	bool m_bRHSAlloc;
	
	Matrix2D<Real[2]> * m_matFourierPsi;
	Matrix2D<Real> * m_matPsi;
	
	void _init();
	
	void _solve(const Real H);
	
	void _copyRHS(const Layer& divergence, Matrix2D<Real>& RHS, const Real dt, const Real rho)
	{
		if(!m_bRHSAlloc)
		{
			RHS=0;
			m_bRHSAlloc = true;
		}
		
		for(int dy=0; dy<divergence.sizeY; dy++)
			for(int dx=0; dx<divergence.sizeX; dx++)
				RHS(dx,dy) = -rho/dt * divergence.read(dx,dy);
	}
	
	void _computeMinusGradP(Layer& minus_grad_pressure)
	{
		const double h= minus_grad_pressure.getH0();
		const double factor = 0.5/h;
		
		for(int dy=0; dy<minus_grad_pressure.sizeY; dy++)
			for(int dx=0; dx<minus_grad_pressure.sizeX; dx++)
			{
				const int dxp = (dx + 1) % minus_grad_pressure.sizeX;
				const int dxm = (dx - 1 + minus_grad_pressure.sizeX) % minus_grad_pressure.sizeX;
				const int dyp = (dy + 1) % minus_grad_pressure.sizeY;
				const int dym = (dy - 1 + minus_grad_pressure.sizeY) % minus_grad_pressure.sizeY;
				
				minus_grad_pressure(dx, dy, 0) = -(m_matPsi->Access(dxp,dy) - m_matPsi->Access(dxm,dy))*factor;
				minus_grad_pressure(dx, dy, 1) = -(m_matPsi->Access(dx,dyp) - m_matPsi->Access(dx,dym))*factor;
			}
	}
	
public:
	PressureSolver_FFTW(const int nX, const int nY) : m_matPsi(NULL), m_matFourierPsi(NULL), m_bReady(false), m_bRHSAlloc(false)
	{
		m_nSize[0] = nX;
		m_nSize[1] = nY;
	}
	
	~PressureSolver_FFTW();
	
	void operator()(const Layer& divergence, Layer& minus_grad_pressure, Real dt, Real rho)
	{
		if (!m_bReady)
		{
			_init();
			m_bReady = true;
		}
		
		_copyRHS(divergence, *m_matPsi, dt, rho);
		_solve(divergence.getH0());
		_computeMinusGradP(minus_grad_pressure);
	}
};

class VelocitySolver_FFTW_BetterAccuracy: public  VelocitySolver_FFTW
{
	Real (* m_fourierU)[2];
	Real (* m_fourierV)[2];
	void * fftw_metadataUV;
	Real * m_destUV[2];
	
	template<typename RealType> 
	void _init(RealType (*&u_fourier)[2], RealType (*&v_fourier)[2], RealType * u_physical, RealType * v_physical);
	
	void _solve(int sizeX, int sizeY, Real H);
public:
	
	VelocitySolver_FFTW_BetterAccuracy(bool bDoublePrecision, Layer& source_vorticity, Layer& dest_velocity);
	
	~VelocitySolver_FFTW_BetterAccuracy();
	
	void operator() (Layer& velocity)
	{
		assert(velocity.sizeX == m_nSize[0]);
		assert(velocity.sizeY == m_nSize[1]);
		
		assert(velocity.getH0() == velocity.getH1());
		_solve(velocity.sizeX, velocity.sizeY, velocity.getH0());
	}
};


class VelocitySolver_Unbounded
{
protected:
	int m_nSize[2], m_nBigSize[2];
	bool m_bReady;
	bool m_bRHSAlloc;
	
	Matrix2D<Real[2]> * m_matFourierGF;
	Matrix2D<Real[2]> * m_matFourierPsi;
	Matrix2D<Real> * m_matPsi;
	
	static inline Real _physGF(int x_, int y_, Real h)
	{
		if (!(x_==0 && y_==0))
		{
			const Real x = (x_+0.0)*h;
			const Real y = (y_+0.0)*h;
			const Real r = sqrt(x*x+y*y);
			return log(r)/(M_PI*2);
		}
		else
		{
			const Real R = h/sqrt(M_PI);
			return (2*log(R)-1)/(4*M_PI);
		}
	}
	
	virtual void _init() =0;
	virtual void _solve()=0;
	
	void _copyRHS(const Layer& vorticity, Matrix2D<Real>& RHS)
	{
		if(!m_bRHSAlloc) 
		{
			RHS=0;
			m_bRHSAlloc = true;
		}

		for(int dy=0; dy<vorticity.sizeY; dy++)
			for(int dx=0; dx<vorticity.sizeX; dx++)
				RHS(dx,dy) = -vorticity.read(dx,dy);
	}
	

	void _computeVelocityFromStreamFunction(Layer& velocity)
	{
		const double h= velocity.getH0();
		const double factor = 0.5/h;
		Matrix2D<Real>& PSI = *m_matPsi;
		
		for(int dy=0; dy<velocity.sizeY; dy++)
			for(int dx=0; dx<velocity.sizeX; dx++)
			{
				/*const int dxp = (dx + 1) % (2*sizeX-2);
				const int dxm = (dx - 1 + 2*sizeX-2) % (2*sizeX-2);
				const int dyp = (dy + 1) % (2*sizeY-2);
				const int dym = (dy - 1 + 2*sizeY-2) % (2*sizeY-2);*/

				const int dxp = (dx + 1) % m_nBigSize[0];
				const int dxm = (dx - 1 + m_nBigSize[0]) % m_nBigSize[0];
				const int dyp = (dy + 1) % m_nBigSize[1];
				const int dym = (dy - 1 + m_nBigSize[1]) % m_nBigSize[1];
				
				velocity(dx, dy, 0) = (PSI(dx,dyp) - PSI(dx,dym))*factor;
				velocity(dx, dy, 1) = -(PSI(dxp,dy) - PSI(dxm,dy))*factor;
			}
	}
	
public:
	VelocitySolver_Unbounded(const int nX, const int nY): 
	m_matFourierGF(NULL), m_matPsi(NULL), m_matFourierPsi(NULL), m_bReady(false), m_bRHSAlloc(false)
	{		
		/*m_nSize[0] = nX;
		m_nSize[1] = nY;
		m_nBigSize[0] = nX*2-2;
		m_nBigSize[1] = nY*2-2;*/
		m_nSize[0] = nX;
		m_nSize[1] = nY;
		m_nBigSize[0] = nX*2;
		m_nBigSize[1] = nY*2;
	}
	
	virtual ~VelocitySolver_Unbounded(){}

	void operator() (const Layer& vorticity, Layer& velocity)
	{
		if (!m_bReady)
		{
			//profiler.getAgent("init").start();
			_init();
			//profiler.getAgent("init").stop();
			m_bReady = true;
		}
		//profiler.getAgent("_copyRHS").start();
		_copyRHS(vorticity, *m_matPsi);
		//profiler.getAgent("_copyRHS").stop();
		_solve();
		//profiler.getAgent("_computeVelFromStreamFunc").start();
		_computeVelocityFromStreamFunction(velocity);
		//profiler.getAgent("_computeVelFromStreamFunc").stop();
		//velocity.MatlabDelCazzo("./matlab/cazzo");
	}
};

class VelocitySolver_Unbounded_FFTW: public VelocitySolver_Unbounded
{
protected:
	
	void _init();
	void _solve();

public:
	
	VelocitySolver_Unbounded_FFTW(const int nX, const int nY): 
		VelocitySolver_Unbounded(nX, nY) {}
	
	~VelocitySolver_Unbounded_FFTW();
};

class PressureSolver_Unbounded_FFTW: public VelocitySolver_Unbounded_FFTW
{
	void _computeMinusGradP(Layer& minus_grad_pressure)
	{
		const double h= minus_grad_pressure.getH0();
		const double factor = 0.5/h;
		//Matrix2D<Real> PSI(this->m_nBigSize[0], this->m_nBigSize[1], (Real *)m_matPsi); // WTF is this?!
		
		/*
		cout << "Checking psi\n";
		for(int dy=0; dy<this->m_nBigSize[1]; dy++)
		{
			for(int dx=0; dx<this->m_nBigSize[0]; dx++)
			{
				//if (isnan(PSI(dx, dy)))
				//	cout << PSI(dx, dy) << " ";
				if (isnan(m_matPsi->Access(dx, dy)))
					cout << m_matPsi->Access(dx, dy) << " ";
			}
		}
		cout << endl;
		cout << h << endl;
		
		//abort();
		//*/
		
		for(int dy=0; dy<minus_grad_pressure.sizeY; dy++)
			for(int dx=0; dx<minus_grad_pressure.sizeX; dx++)
			{
				const int dxp = (dx + 1) % (2*minus_grad_pressure.sizeX-2);
				const int dxm = (dx - 1 + 2*minus_grad_pressure.sizeX-2) % (2*minus_grad_pressure.sizeX-2);
				const int dyp = (dy + 1) % (2*minus_grad_pressure.sizeY-2);
				const int dym = (dy - 1 + 2*minus_grad_pressure.sizeY-2) % (2*minus_grad_pressure.sizeY-2);
				
				minus_grad_pressure(dx, dy, 0) = (m_matPsi->Access(dxp,dy) - m_matPsi->Access(dxm,dy))*factor;
				minus_grad_pressure(dx, dy, 1) = (m_matPsi->Access(dx,dyp) - m_matPsi->Access(dx,dym))*factor;
				
				//if (isnan(minus_grad_pressure(dx, dy, 0)))
				//	cout << minus_grad_pressure(dx, dy, 0) << "\t" << PSI(dxp,dy) << "\t" << PSI(dxm,dy) << "\t" << dxp << "\t" << dxm << "\t" << dy << "\t" << this->m_nBigSize[0] << "\t" << this->m_nBigSize[1] << "\t" << factor << endl;
			}
		/*
		for(int dy=0; dy<minus_grad_pressure.sizeY; dy++)
		{
			for(int dx=0; dx<minus_grad_pressure.sizeX; dx++)
			{
				cout << minus_grad_pressure(dx, dy, 0) << " ";
			}
			cout << endl;
		}
		cout << endl;
		*/
		//abort();
	}	
public:
	PressureSolver_Unbounded_FFTW(const int nX, const int nY):
		VelocitySolver_Unbounded_FFTW(nX, nY)
	{
	}
	
	void operator() (const Layer& divergence, Layer& minus_grad_pressure)
	{
		if (!m_bReady)
		{
			_init();
			m_bReady = true;
		}
		
		_copyRHS(divergence, *m_matPsi);
		_solve();
		_computeMinusGradP(minus_grad_pressure);
	}
};

class VelocitySolver_Unbounded_FFTW_BetterAccuracy: public VelocitySolver_Unbounded_FFTW
{
	//_FORBID_COPIES(VelocitySolver_Unbounded_FFTW_BetterAccuracy);
	
	Matrix2D<Real> * m_matVelocity[2];
	Matrix2D<Real[2]> * m_matFourierVelocity[2];
	
	void _init();
	void _solve();
	
	template <int iDim>
	void _copyVelocity(Layer& velocity)
	{
		Matrix2D<Real>& VEL = *m_matVelocity[iDim];
		for(int dy=0; dy<velocity.sizeY; dy++)
			for(int dx=0; dx<velocity.sizeX; dx++)
				velocity(dx, dy, iDim) = VEL(dx,dy);
	}
	
public:
	VelocitySolver_Unbounded_FFTW_BetterAccuracy(const int nX, const int nY): 
		VelocitySolver_Unbounded_FFTW(nX, nY)
	{
		m_matVelocity[0] = NULL;
		m_matVelocity[1] = NULL;
		m_matFourierVelocity[0] = NULL;
		m_matFourierVelocity[1] = NULL;
		
		VelocitySolver_Unbounded_FFTW::_init();
		_init();
	}
	
	~VelocitySolver_Unbounded_FFTW_BetterAccuracy();
	
	void operator() (const Layer& vorticity, Layer& velocity)
	{
		//vorticity.MatlabDelCazzo("./matlab/input");
		_copyRHS(vorticity, *m_matVelocity[0]);
		_solve();
		_copyVelocity<0>(velocity);
		_copyVelocity<1>(velocity);
		//velocity.MatlabDelCazzo("./matlab/cazzo");
	}		
};

