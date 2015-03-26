/*
 *  Operators_DFT.cpp
 *  VM2D
 *
 *  Created by Diego Rossinelli on 2/10/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include <fftw3.h>
#include <math.h>
#include "Profiler.h"

extern const bool bProfile;
extern Profiler profiler;

#include "Operators_DFT.h"
#include "InterfaceFFTW.h"

struct MetadataFFWT_float
{
	fftwf_plan plan_forward;
	fftwf_plan plan_backward;
	fftwf_plan velocityplan_backward[2];
	
};
struct MetadataFFWT_double
{
	fftw_plan plan_forward;
	fftw_plan plan_backward;
	fftw_plan velocityplan_backward[2];
};

void VelocitySolver_FFTW::_dispose()
{
	if (m_bDoublePrecision)
	{
		fftw_free((void *)data_fourier);
		fftw_destroy_plan(((MetadataFFWT_double *)fftw_metadata)->plan_forward);
		fftw_destroy_plan(((MetadataFFWT_double *)fftw_metadata)->plan_backward);
	}
	else
	{
		fftwf_free((void *)data_fourier);
		fftwf_destroy_plan(((MetadataFFWT_float *)fftw_metadata)->plan_forward);
		fftwf_destroy_plan(((MetadataFFWT_float *)fftw_metadata)->plan_backward);
	}
}

template<>
void VelocitySolver_FFTW::_init<double>(double * vorticity, double (*&data_fourier)[2], double * psi_physical)
{
	MetadataFFWT_double * metadata = new MetadataFFWT_double;
	fftw_metadata = metadata;
	
	//allocate data
	double  * datain = (double*)fftw_malloc(2*sizeof(double)*(m_nSize[0]/2+1)*(m_nSize[1]));
	data_fourier = (double (*)[2])datain;
	
	//create plans
	metadata->plan_forward = fftw_plan_dft_r2c_2d(m_nSize[1], m_nSize[0], vorticity, data_fourier, FFTW_MEASURE);
	metadata->plan_backward = fftw_plan_dft_c2r_2d(m_nSize[1], m_nSize[0], data_fourier, psi_physical, FFTW_MEASURE);
}

template<>
void VelocitySolver_FFTW::_init<float>(float * vorticity, float (*&data_fourier)[2], float * psi_physical)
{
	MetadataFFWT_float * metadata = new MetadataFFWT_float;
	fftw_metadata = metadata;
	
	//allocate data
	float  * datain = (float*)fftwf_malloc(2*sizeof(float)*(m_nSize[0]/2+1)*(m_nSize[1]));
	data_fourier = (float (*)[2])datain;
	
	//create plans
	metadata->plan_forward = fftwf_plan_dft_r2c_2d(m_nSize[1], m_nSize[0], vorticity, data_fourier, FFTW_MEASURE);
	metadata->plan_backward = fftwf_plan_dft_c2r_2d(m_nSize[1], m_nSize[0], data_fourier, psi_physical, FFTW_MEASURE);
}

VelocitySolver_FFTW::VelocitySolver_FFTW(bool bDoublePrecision, Layer& vorticity):
psi_physical(NULL), fftw_metadata(NULL), m_bDoublePrecision(bDoublePrecision)
{
	m_nSize[0] = vorticity.sizeX;
	m_nSize[1] = vorticity.sizeY;
	
	Layer * psi_layer = new Layer(vorticity.sizeX, vorticity.sizeY, 1);
	psi_physical = (Real *)psi_layer->data;
	
	_init((Real*)vorticity.data, data_fourier, psi_physical);
}

void VelocitySolver_FFTW::_solve_phi_core(Real (*&data_fourier)[2], Real H)
{
	const int nX = m_nSize[0]/2+1;
	const int nY = m_nSize[1];
	
	const Real H2 = (H*H);
	const Real invNX = 1.0/m_nSize[0];
	const Real invNY = 1.0/m_nSize[1];
	
	data_fourier[0][0] = 0;
	data_fourier[0][1] = 0;
	for(int iy=0; iy<nY; iy++)
		for(int ix=0; ix<nX; ix++)
		{
			const Real denom = 2*cos(2*M_PI*ix*invNX) -4 +2*cos(2*M_PI*iy*invNY);
			const Real laplace_factor = (ix==0 && iy==0)? 0 : (-H2/denom);
			
			data_fourier[nX*iy + ix][0] *= laplace_factor;
			data_fourier[nX*iy + ix][1] *= laplace_factor;
		}
}

void VelocitySolver_FFTW::_solve_phi(int sizeX, int sizeY, Real H)
{
	if (m_bDoublePrecision)
		fftw_execute(((MetadataFFWT_double *)fftw_metadata)->plan_forward);
	else
		fftwf_execute(((MetadataFFWT_float *)fftw_metadata)->plan_forward);
	
	_solve_phi_core(data_fourier, H);
	
	if (m_bDoublePrecision)
		fftw_execute(((MetadataFFWT_double *)fftw_metadata)->plan_backward);
	else
		fftwf_execute(((MetadataFFWT_float *)fftw_metadata)->plan_backward);
	
	
	const Real invNXNY = 1.0/(sizeY*sizeX);
	for(int dy=0; dy<sizeY; dy++)
		for(int dx=0; dx<sizeX; dx++)
			psi_physical[dx +  sizeX*dy] *= invNXNY;
}


void PressureSolver_FFTW::_init()
{
	InterfaceFFTW::allocateMatrix(m_matFourierPsi, m_nSize[0]/2+1, m_nSize[1]);
	InterfaceFFTW::allocateMatrix(m_matPsi, m_nSize[0], m_nSize[1]);
	
	InterfaceFFTW::createPlanForward("RHS: Real->Fourier", *m_matPsi, *m_matFourierPsi);
	InterfaceFFTW::createPlanBackward("Psi: Fourier->Real", *m_matFourierPsi, *m_matPsi);
}

void PressureSolver_FFTW::_solve(const Real H)
{
	const Real invNX = 1.0/m_nSize[0];
	const Real invNY = 1.0/m_nSize[1];
	const Real H2 = (H*H);
	const int nX = m_nSize[0]/2+1;
	const int nY = m_nSize[1];
	
	// the 2 components are the real and imaginary parts
	Matrix2D<Real[2]>& PSI = *m_matFourierPsi;
	
	//profiler.getAgent("FFTW-F").start();
	InterfaceFFTW::executePlan("RHS: Real->Fourier");
	//profiler.getAgent("FFTW-F").stop();
	
	//profiler.getAgent("FFTW-convolution").start();
	const Real WcoeffX[2] = {cos(2*M_PI/nX),sin(2*M_PI/nX)};
	const Real WcoeffY[2] = {cos(2*M_PI/nY),sin(2*M_PI/nY)};
	Real Wx[2] = {1.,1.};
	Real Wy[2] = {1.,1.};
	for(int iy=0; iy<nY; iy++)
	{
		for(int ix=0; ix<nX; ix++)
		{
			Real denom[2] = {4.,0.};
			Real denomWx = Wx[0]*Wx[0] + Wx[1]*Wx[1];
			Real denomWy = Wx[0]*Wx[0] + Wx[1]*Wx[1];
			denom[0] -= Wx[0] + Wx[0]/denomWx + Wy[0] + Wy[0]/denomWy;
			denom[1] -= Wx[1] - Wx[1]/denomWx + Wy[1] - Wy[1]/denomWy;
			
			if (denom[0]!=0 && denom[1]!=0)
			{
				Real denomLength = denom[0]*denom[0]+denom[1]*denom[1];
				Real invDenom[2] = {denom[0]/denomLength,-denom[1]/denomLength};
				
				PSI(ix,iy)[0] *= H*H * invDenom[0];
				PSI(ix,iy)[1] *= H*H * invDenom[1];
			}
			
			Wx[0] = Wx[0]*WcoeffX[0] - Wx[1]*WcoeffX[1];
			Wx[1] = Wx[0]*WcoeffX[1] - Wx[1]*WcoeffX[0];
		}
		Wy[0] = Wy[0]*WcoeffY[0] - Wy[1]*WcoeffY[1];
		Wy[1] = Wy[0]*WcoeffY[1] - Wy[1]*WcoeffY[0];
	}
	//profiler.getAgent("FFTW-convolution").stop();
	
	//profiler.getAgent("FFTW-B").start();
	InterfaceFFTW::executePlan("Psi: Fourier->Real");
	//profiler.getAgent("FFTW-B").stop();
}

PressureSolver_FFTW::~PressureSolver_FFTW()
{
	InterfaceFFTW::erasePlan("RHS: Real->Fourier");
	InterfaceFFTW::erasePlan("Psi: Fourier->Real");
	InterfaceFFTW::deallocateMatrix(m_matFourierPsi);
	InterfaceFFTW::deallocateMatrix(m_matPsi);
}

template <>
void VelocitySolver_FFTW_BetterAccuracy::_init<double>(double (*&u_fourier)[2], double (*&v_fourier)[2], double * u_physical, double * v_physical)
{
	vector<fftw_plan> & v = *(new vector<fftw_plan>(2));
	fftw_metadataUV = &v;
	
	u_fourier = (double (*)[2])fftw_malloc(2*sizeof(double)*(m_nSize[0]/2+1)*(m_nSize[1]));
	v_fourier = (double (*)[2])fftw_malloc(2*sizeof(double)*(m_nSize[0]/2+1)*(m_nSize[1]));
	
	v[0] = fftw_plan_dft_c2r_2d(m_nSize[1], m_nSize[0], u_fourier, u_physical, FFTW_MEASURE);
	v[1] = fftw_plan_dft_c2r_2d(m_nSize[1], m_nSize[0], v_fourier, v_physical, FFTW_MEASURE);
}


template <>
void VelocitySolver_FFTW_BetterAccuracy::_init<float>(float (*&u_fourier)[2], float (*&v_fourier)[2], float * u_physical, float * v_physical)
{
	vector<fftwf_plan> & v = *(new vector<fftwf_plan>(2));
	fftw_metadataUV = &v;
	
	u_fourier = (float (*)[2])fftwf_malloc(2*sizeof(float)*(m_nSize[0]/2+1)*(m_nSize[1]));
	v_fourier = (float (*)[2])fftwf_malloc(2*sizeof(float)*(m_nSize[0]/2+1)*(m_nSize[1]));
	
	v[0] = fftwf_plan_dft_c2r_2d(m_nSize[1], m_nSize[0], u_fourier, u_physical, FFTW_MEASURE);
	v[1] = fftwf_plan_dft_c2r_2d(m_nSize[1], m_nSize[0], v_fourier, v_physical, FFTW_MEASURE);
}

VelocitySolver_FFTW_BetterAccuracy::VelocitySolver_FFTW_BetterAccuracy(bool bDoublePrecision, Layer& source_vorticity, Layer& dest_velocity):
VelocitySolver_FFTW(bDoublePrecision, source_vorticity)
{
	Real * ptrU = dest_velocity.template getPlane<0>();
	Real * ptrV = dest_velocity.template getPlane<1>();
	
	assert(((unsigned int)(ptrU) & 0x0F) == 0);
	assert(((unsigned int)(ptrV) & 0x0F) == 0);
	
	m_destUV[0] = ptrU;
	m_destUV[1] = ptrV;
	
	_init(m_fourierU, m_fourierV, m_destUV[0], m_destUV[1]);
}


VelocitySolver_FFTW_BetterAccuracy::~VelocitySolver_FFTW_BetterAccuracy()
{
	if (m_bDoublePrecision)
	{
		fftw_free((void *)m_fourierU);
		fftw_free((void *)m_fourierV);
		vector<fftw_plan> & v = *((vector<fftw_plan>*)fftw_metadataUV);
		fftw_destroy_plan(v[0]);
		fftw_destroy_plan(v[1]);
	}
	else
	{
		fftwf_free((void *)m_fourierU);
		fftwf_free((void *)m_fourierV);
		vector<fftwf_plan> & v = *((vector<fftwf_plan>*)fftw_metadataUV);
		fftwf_destroy_plan(v[0]);
		fftwf_destroy_plan(v[1]);
	}
}

void VelocitySolver_FFTW_BetterAccuracy::_solve(int sizeX, int sizeY, Real H)
{
	if (m_bDoublePrecision)
		fftw_execute(((MetadataFFWT_double *)fftw_metadata)->plan_forward);
	else
		fftwf_execute(((MetadataFFWT_float *)fftw_metadata)->plan_forward);
	
	_solve_phi_core(data_fourier, H);
	
	const int nX = m_nSize[0]/2+1;
	const int nY = m_nSize[1];
	const Real pi2 = 2*M_PI;
	
	for(int dy=0; dy<nY; dy++)
		for(int dx=0; dx<nX; dx++)
		{
			const int ky =  ((dy + nY/2) % nY) - nY/2;
			m_fourierU[dx +  nX*dy][0] = -data_fourier[dx +  nX*dy][1]*ky*pi2;
			m_fourierU[dx +  nX*dy][1] = data_fourier[dx +  nX*dy][0]*ky*pi2;
		}
	
	for(int dy=0; dy<nY; dy++)
		for(int dx=0; dx<nX; dx++)
		{
			const int kx = dx;
			m_fourierV[dx +  nX*dy][0] = data_fourier[dx +  nX*dy][1]*kx*pi2;
			m_fourierV[dx +  nX*dy][1] = -data_fourier[dx +  nX*dy][0]*kx*pi2;
		}
	
	
	if (m_bDoublePrecision)
	{
		vector<fftw_plan> & v = *((vector<fftw_plan>*)fftw_metadataUV);
		fftw_execute(v[0]);
		fftw_execute(v[1]);
	}
	else
	{
		vector<fftwf_plan> & v = *((vector<fftwf_plan>*)fftw_metadataUV);
		fftwf_execute(v[0]);
		fftwf_execute(v[1]);
	}
	
	const Real invNXNY = 1./(sizeY*sizeX);
	
	Real * u = m_destUV[0];
	for(int dy=0; dy<sizeY; dy++)
		for(int dx=0; dx<sizeX; dx++)
			u[dx +  sizeX*dy] *= invNXNY;
	
	Real * v = m_destUV[1];
	for(int dy=0; dy<sizeY; dy++)
		for(int dx=0; dx<sizeX; dx++)
			v[dx +  sizeX*dy] *= invNXNY;
}


VelocitySolver_Unbounded_FFTW::~VelocitySolver_Unbounded_FFTW()
{
	InterfaceFFTW::erasePlan("RHS: Real->Fourier");
	InterfaceFFTW::erasePlan("Psi: Fourier->Real");
	InterfaceFFTW::deallocateMatrix(m_matFourierGF);
	InterfaceFFTW::deallocateMatrix(m_matFourierPsi);
	InterfaceFFTW::deallocateMatrix(m_matPsi);
}

void VelocitySolver_Unbounded_FFTW::_init()
{
	InterfaceFFTW::allocateMatrix(m_matFourierGF, m_nBigSize[0]/2+1, m_nBigSize[1]);
	InterfaceFFTW::allocateMatrix(m_matFourierPsi, m_nBigSize[0]/2+1, m_nBigSize[1]);
	InterfaceFFTW::allocateMatrix(m_matPsi, m_nBigSize[0], m_nBigSize[1]);
	
	InterfaceFFTW::createPlanForward("RHS: Real->Fourier", *m_matPsi, *m_matFourierPsi);
	InterfaceFFTW::createPlanBackward("Psi: Fourier->Real", *m_matFourierPsi, *m_matPsi);
	
	{
		Matrix2D<Real> * matGF  = NULL;
		InterfaceFFTW::allocateMatrix(matGF, m_nBigSize[0], m_nBigSize[1]);
		InterfaceFFTW::createPlanForward("GF: Real->Fourier", *matGF, *m_matFourierGF);
		
		const Real H = 1./m_nSize[0];
		for(int dy=0; dy<m_nBigSize[1]; dy++)
			for(int dx=0; dx<m_nBigSize[0]; dx++)
				matGF->Access(dx, dy) = _physGF(min(dx, m_nBigSize[0]-dx), min(dy, m_nBigSize[1]-dy), H);
		
		InterfaceFFTW::executePlan("GF: Real->Fourier");
		InterfaceFFTW::erasePlan("GF: Real->Fourier");
		InterfaceFFTW::deallocateMatrix(matGF);
	}
}

void VelocitySolver_Unbounded_FFTW::_solve()
{
	const Real h = 1./m_nSize[0];
	const int nX = m_nBigSize[0]/2+1;
	const int nY = m_nBigSize[1];
	
	Matrix2D<Real[2]>& PSI = *m_matFourierPsi;
	Matrix2D<Real[2]>& GF = *m_matFourierGF;
	
	//profiler.getAgent("FFTW-F").start();
	InterfaceFFTW::executePlan("RHS: Real->Fourier");
	//profiler.getAgent("FFTW-F").stop();
	
	//profiler.getAgent("FFTW-convolution").start();
	for(int dy=0; dy<nY; dy++)
		for(int dx=0; dx<nX; dx++)
		{
			const Real a[2] = {PSI(dx, dy)[0], PSI(dx, dy)[1] };
			const Real b[2] = {GF(dx, dy)[0], GF(dx, dy)[1] };
			
			PSI(dx, dy)[0] = (a[0]*b[0]-a[1]*b[1]);
			PSI(dx, dy)[1] = (a[1]*b[0]+a[0]*b[1]);
		}
	//profiler.getAgent("FFTW-convolution").stop();
	
	//profiler.getAgent("FFTW-B").start();
	InterfaceFFTW::executePlan("Psi: Fourier->Real");
	(*m_matPsi) *= h*h/(m_nBigSize[0]*m_nBigSize[1]);
	//profiler.getAgent("FFTW-B").stop();
	
	
	
	/*
	for(int dy=0; dy<m_nBigSize[1]; dy++)
	{
		for(int dx=0; dx<m_nBigSize[0]; dx++)
		{
			if (isnan(m_matPsi->Access(dx, dy)))
				cout << m_matPsi->Access(dx, dy) << " ";
		}
	}
	cout << endl;
	
	cout << "Done with fourier solver\n";
	//abort();
	//*/
}

void VelocitySolver_Unbounded_FFTW_BetterAccuracy::_init()
{
	InterfaceFFTW::deallocateMatrix(m_matPsi);
	
	InterfaceFFTW::allocateMatrix(m_matFourierVelocity[0], m_nBigSize[0]/2+1, m_nBigSize[1]);
	InterfaceFFTW::allocateMatrix(m_matFourierVelocity[1], m_nBigSize[0]/2+1, m_nBigSize[1]);
	InterfaceFFTW::allocateMatrix(m_matVelocity[0], m_nBigSize[0], m_nBigSize[1]);
	InterfaceFFTW::allocateMatrix(m_matVelocity[1], m_nBigSize[0], m_nBigSize[1]);
	
	InterfaceFFTW::createPlanForward("tricky RHS: Real->Fourier", *m_matVelocity[0], *m_matFourierPsi);
	InterfaceFFTW::createPlanBackward("U: Fourier->Real", *m_matFourierVelocity[0], *m_matVelocity[0]);
	InterfaceFFTW::createPlanBackward("V: Fourier->Real", *m_matFourierVelocity[1], *m_matVelocity[1]);
}

void VelocitySolver_Unbounded_FFTW_BetterAccuracy::_solve()
{
	const Real h = 1./m_nSize[0];
	
	const double factorU = 2*M_PI*h/m_nBigSize[1];
	const double factorV = 2*M_PI*h/m_nBigSize[0];
	
	const int nX = m_nBigSize[0]/2+1;
	const int nY = m_nBigSize[1];
	
	Matrix2D<Real[2]>& PSI = *m_matFourierPsi;
	Matrix2D<Real[2]>& GF = *m_matFourierGF;
	Matrix2D<Real[2]>& U = *m_matFourierVelocity[0];
	Matrix2D<Real[2]>& V = *m_matFourierVelocity[1];
	
	InterfaceFFTW::executePlan("tricky RHS: Real->Fourier");
	
	for(int dy=0; dy<nY; dy++)
		for(int dx=0; dx<nX; dx++)
		{
			const Real a[2] = {PSI(dx, dy)[0], PSI(dx, dy)[1] };
			const Real b[2] = {GF(dx, dy)[0], GF(dx, dy)[1] };
			
			PSI(dx, dy)[0] = (a[0]*b[0]-a[1]*b[1]);
			PSI(dx, dy)[1] = (a[1]*b[0]+a[0]*b[1]);
		}
	
	for(int dy=0; dy<nY; dy++)
		for(int dx=0; dx<nX; dx++)
		{
			const int ky =  ((dy + nY/2) % nY) - nY/2;
			
			U(dx, dy)[0] = -ky*factorU*PSI(dx, dy)[1];
			U(dx, dy)[1] = ky*factorU*PSI(dx, dy)[0];
		}
	
	for(int dy=0; dy<nY; dy++)
		for(int dx=0; dx<nX; dx++)
		{
			const int kx = dx;
			
			V(dx, dy)[0] = kx*factorV*PSI(dx, dy)[1];
			V(dx, dy)[1] = -kx*factorV*PSI(dx, dy)[0];
		}
	
	InterfaceFFTW::executePlan("U: Fourier->Real");
	InterfaceFFTW::executePlan("V: Fourier->Real");
	
	(*m_matVelocity[0]) *= 1./(m_nBigSize[0]*m_nBigSize[1]);
	(*m_matVelocity[1]) *= 1./(m_nBigSize[0]*m_nBigSize[1]);
}

VelocitySolver_Unbounded_FFTW_BetterAccuracy::~VelocitySolver_Unbounded_FFTW_BetterAccuracy()
{
	InterfaceFFTW::deallocateMatrix(m_matFourierVelocity[0]);
	InterfaceFFTW::deallocateMatrix(m_matFourierVelocity[1]);
	InterfaceFFTW::deallocateMatrix(m_matVelocity[0]);
	InterfaceFFTW::deallocateMatrix(m_matVelocity[1]);
	
	InterfaceFFTW::erasePlan("tricky RHS: Real->Fourier");
	InterfaceFFTW::erasePlan("U: Fourier->Real");
	InterfaceFFTW::erasePlan("V: Fourier->Real");
} 

