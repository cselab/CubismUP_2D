/*
 *  InterfaceFFTW.h
 *  VM2D
 *
 *  Created by Diego Rossinelli on 2/25/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include <fftw3.h>

#include <map>
using namespace std;

class InterfaceFFTW
{
	static const int bMultithreading = true;
	static const int nThreads = NTHREADS;
	
	static InterfaceFFTW singletone;
	static map<string, bool > mapDoublePrecision;
	static map<string, void *> mapPlans;
	
	template<typename T> static void _deallocateRegion(T* t);
	template<typename T> static void _allocateRegion(int nX, int nY, T*& t);
	template<typename T> static void * _createPlan(T* real_data, T(*complex_data)[2], const int nX, const int nY, int sign);
	
	InterfaceFFTW()
	{
		fftwf_init_threads();
		fftwf_plan_with_nthreads(nThreads);
		printf("FFTW initialized with %d threads\n", nThreads);
	}
public:
	template<typename T>
	static void allocateMatrix(Matrix2D<T>*& mat, int nX, int nY)
	{
		assert(mat==NULL);
		
		T* t = NULL;
		_allocateRegion(nX, nY, t);
		assert(t!=NULL);
		
		mat = new Matrix2D<T>(nX, nY, t);
		assert(mat !=NULL);
	}

	template<typename T>
	static void deallocateMatrix(Matrix2D<T>*& mat)
	{
		if (mat == NULL) return;
		//assert(mat!=NULL);
		
		T* t = &mat->LinAccess(0);
		_deallocateRegion(t);
		
		delete mat; 
		mat = NULL;
	}
	
	template<typename T>
	static void createPlanForward(string sPlanName, Matrix2D<T>& matIn, Matrix2D<T[2]>& matOut)
	{
		assert(sizeof(T)==8 || sizeof(T)==4);
		assert(matIn.getSize()[0]/2+1 == matOut.getSize()[0]);
		assert(matIn.getSize()[1] == matOut.getSize()[1]);
		
		const map<string, void *>::const_iterator it = mapPlans.find(sPlanName);
		assert(it == mapPlans.end());
		
		T* in = &matIn.LinAccess(0);
		T(* out)[2] = &matOut.LinAccess(0);
		
		mapPlans[sPlanName] = _createPlan(in,out, matIn.getSize()[0], matIn.getSize()[1],+1);
		mapDoublePrecision[sPlanName] = (sizeof(T)==8);
	}
	
	template<typename T>
	static void createPlanBackward(string sPlanName, Matrix2D<T[2]>& matIn, Matrix2D<T>& matOut)
	{
		assert(matIn.getSize()[0] == matOut.getSize()[0]/2+1);
		assert(matIn.getSize()[1] == matOut.getSize()[1]);
		
		T(*in)[2] = &matIn.LinAccess(0);
		T* out = &matOut.LinAccess(0);
		mapPlans[sPlanName] = _createPlan(out, in, matOut.getSize()[0], matOut.getSize()[1],-1);
	}
	

	static void erasePlan(string sPlanName)
	{
		const map<string, void *>::const_iterator it = mapPlans.find(sPlanName);
		
		assert(it != mapPlans.end());
		
		const bool bDoublePrecision = mapDoublePrecision[sPlanName];
		
		if(bDoublePrecision)
			fftw_destroy_plan(*(fftw_plan *)it->second);
		else
			fftwf_destroy_plan(*(fftwf_plan *)it->second);
		
		delete (char *)it->second;
		
		
		mapDoublePrecision.erase(sPlanName);
		mapPlans.erase(sPlanName);
	}
	
	static void executePlan(string sPlanName) 
	{
		const map<string, void *>::const_iterator it = mapPlans.find(sPlanName);
		
		assert(it != mapPlans.end());
		
		const bool bDoublePrecision = mapDoublePrecision[sPlanName];
		
		if(bDoublePrecision)
			fftw_execute(*(fftw_plan *)it->second);
		else
			fftwf_execute(*(fftwf_plan *)it->second);
	}
};

InterfaceFFTW InterfaceFFTW::singletone;
map<string, bool > InterfaceFFTW::mapDoublePrecision;
map<string, void *> InterfaceFFTW::mapPlans;

template <>
void InterfaceFFTW::_allocateRegion(int nX, int nY, double*& t)
{
	t = (double*)fftw_malloc(sizeof(double)*nX*nY);
}

template <>
void InterfaceFFTW::_allocateRegion(int nX, int nY, float*& t)
{
	t = (float*)fftwf_malloc(sizeof(float)*nX*nY);
}

template <>
void InterfaceFFTW::_allocateRegion(int nX, int nY, double (*& t)[2])
{
	t = (double (*)[2])fftw_malloc(sizeof(double[2])*nX*nY);
}

template <>
void InterfaceFFTW::_allocateRegion(int nX, int nY, float (*& t)[2])
{
	t = (float (*)[2])fftwf_malloc(sizeof(float[2])*nX*nY);
}

template <> void InterfaceFFTW::_deallocateRegion(float * t) {fftwf_free(t);}
template <> void InterfaceFFTW::_deallocateRegion(double * t) {fftw_free(t);}
template <> void InterfaceFFTW::_deallocateRegion(float (* t)[2]) {fftwf_free((float*)t);}
template <> void InterfaceFFTW::_deallocateRegion(double (* t)[2]) {fftw_free((double*)t);}

template<> void * InterfaceFFTW::_createPlan<float>(float* real_data, float (*complex_data)[2], const int nX, const int nY, int sign)
{
	fftwf_plan * plan = new fftwf_plan;
	
	if (sign==1)
		*plan = fftwf_plan_dft_r2c_2d(nY, nX, real_data, complex_data, FFTW_MEASURE);
	else if (sign==-1)
		*plan = fftwf_plan_dft_c2r_2d(nY, nX, complex_data, real_data, FFTW_MEASURE);
	else abort();
	
	return plan;
}

template<> void * InterfaceFFTW::_createPlan<double>(double* real_data, double (*complex_data)[2], const int nX, const int nY, int sign)
{
	fftw_plan * plan = new fftw_plan;
	
	if (sign==1)
		*plan = fftw_plan_dft_r2c_2d(nY, nX, real_data, complex_data, FFTW_MEASURE);
	else if (sign==-1)
		*plan = fftw_plan_dft_c2r_2d(nY, nX, complex_data, real_data, FFTW_MEASURE);
	else abort();
	
	return plan;
}

