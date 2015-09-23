#pragma once

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>

#include <fftw3.h>

#include "common.h"

#ifndef _SP_COMP_
// double
typedef fftw_complex mycomplex;
typedef fftw_plan myplan;
#else // _SP_COMP_
// float
typedef fftwf_complex mycomplex;
typedef fftwf_plan myplan;
#endif // _SP_COMP_

using namespace std;

#include <BlockInfo.h>
#include <Profiler.h>

class FFTWBase
{
	static int registered_objects;
	static bool initialized; //a la singleton
	
	static void _setup(const int desired_threads)
	{
		if (!initialized)
		{
			initialized = true;
			
#ifndef _SP_COMP_
			const int retval = fftw_init_threads();
#else // _SP_COMP_
			const int retval = fftwf_init_threads();
#endif // _SP_COMP_
			if(retval==0)
			{
				cout << "FFTWBase::setup(): Oops the call to fftw_init_threads() returned zero. Aborting\n";
				abort();
			}
			else
#ifndef _SP_COMP_
				fftw_plan_with_nthreads(desired_threads);
#else // _SP_COMP_
			fftwf_plan_with_nthreads(desired_threads);
#endif // _SP_COMP_
		}
		
		registered_objects++;
	}
	
public:
	
	FFTWBase(const int desired_threads) { _setup(desired_threads); }
	
	static void dispose()
	{
		registered_objects--;
		
		if (registered_objects == 0)
		{
#ifndef _SP_COMP_
			fftw_cleanup_threads();
#else // _SP_COMP_
			fftwf_cleanup_threads();
#endif // _SP_COMP_
		}
	}
};

template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW : FFTWBase
{
	Profiler profiler;
	
protected:
	typedef typename TGrid::BlockType BlockType;
	
	bool initialized;
	
	myplan fwd, bwd;
	Real * data; // rhs in _setup, out in cub2fftw and fftw2cub
	
protected:
	
	virtual void _setup(Real *& rhs, const size_t nx, const size_t ny)
	{
		if (!initialized)
		{
			initialized = true;
#ifndef _SP_COMP_
			rhs = fftw_alloc_real(2*nx*(ny/2+1));
			
			fwd = fftw_plan_dft_r2c_2d(nx, ny, rhs, (mycomplex *)rhs, FFTW_MEASURE);
			bwd = fftw_plan_dft_c2r_2d(nx, ny, (mycomplex *)rhs, rhs, FFTW_MEASURE);
#else // _SP_COMP_
			rhs = fftwf_alloc_real(2*nx*(ny/2+1));
			
			fwd = fftwf_plan_dft_r2c_2d(nx, ny, rhs, (mycomplex *)rhs, FFTW_MEASURE);
			bwd = fftwf_plan_dft_c2r_2d(nx, ny, (mycomplex *)rhs, rhs, FFTW_MEASURE);
#endif // _SP_COMP_
		}
	}
	
	void _cub2fftw(TGrid& grid, Real * out, const size_t nx, const size_t ny, const size_t ny_hat) const
	{
		vector<BlockInfo> infos = grid.getBlocksInfo();
		
		const size_t N = infos.size();
		
#pragma omp parallel for
		for(int i=0; i<N; ++i)
		{
			const BlockInfo info = infos[i];
			BlockType& b = *(BlockType*)infos[i].ptrBlock;
			
			// indexing for fftw - [x][y] (y is the fastest running index), block part of the index
			const size_t offset = BlockType::sizeY*info.index[1] + (ny/2+1)*2*BlockType::sizeX*info.index[0];
			
			for(int iy=0; iy<BlockType::sizeY; iy++)
				for(int ix=0; ix<BlockType::sizeX; ix++)
				{
					// indexing for fftw - [x][y] (y is the fastest running index), local part of the index
					const size_t dest_index = TStreamer::channels*(offset + iy + 2*(ny/2+1)*ix);
					
					assert(dest_index>=0 && dest_index<nx*ny_hat*2); // linking error with openmp: http://lists.apple.com/archives/xcode-users/2011/Oct/msg00252.html
					TStreamer::operate(b(ix,iy), &out[dest_index]);
				}
		}
	}
	
	virtual void _solve(mycomplex * in_out, const size_t nx, const size_t ny, const size_t ny_hat, const Real norm_factor, const Real h)
	{
		if (TStreamer::channels != 1)
		{
			cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
			abort();
		}
		
		const Real h2 = h*h;
		const Real factor = h2*norm_factor;
		
#pragma omp parallel for
		for(int i=0; i<nx; ++i)
			for(int j = 0; j<ny_hat; ++j)
			{
				const int linidx = i*ny_hat+j;
				assert(linidx >=0 && linidx<nx*ny_hat); // linking error with openmp
				
				// based on the 5 point stencil in 1D (h^4 error)
				const Real denom = 32.*(cos(2.*M_PI*i/nx) + cos(2.*M_PI*j/ny)) - 2.*(cos(4.*M_PI*i/nx) + cos(4.*M_PI*j/ny)) - 60.;
				const Real inv_denom = (denom==0)? 0.:1./denom;
				const Real fatfactor = 12. * inv_denom * factor;
				
				// based on the 3 point stencil in 1D (h^2 error)
				//const Real denom = 2.*(cos(2.*M_PI*i/nx) + cos(2.*M_PI*j/ny)) - 4.;
				//const Real inv_denom = (denom==0)? 0.:1./denom;
				//const Real fatfactor = inv_denom * factor;
				
				// this is to check the transform only
				//const Real fatfactor = norm_factor;
				
				in_out[linidx][0] *= fatfactor;
				in_out[linidx][1] *= fatfactor;
			}
		
		//this is sparta!
		in_out[0][0] = in_out[0][1] = 0;
	}
	
	virtual void _solveSpectral(mycomplex * in_out, const size_t nx, const size_t ny, const size_t ny_hat, const Real norm_factor, const Real h)
	{
		if (TStreamer::channels != 1)
		{
			cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
			abort();
		}
		
		const Real waveFactX = 2.0*M_PI/(nx*h);
		const Real waveFactY = 2.0*M_PI/(ny*h);
#pragma omp parallel for
		for(int i=0; i<nx; ++i)
			for(int j = 0; j<ny_hat; ++j)
			{
				const int linidx = i*ny_hat+j;
				assert(linidx >=0 && linidx<nx*ny_hat); // linking error with openmp
				
				// wave number
				const int kx = (i <= nx/2) ? i : -(nx-i);
				const int ky = (j <= ny/2) ? j : -(ny-j);
				const Real rkx = kx*waveFactX;
				const Real rky = ky*waveFactY;
				const Real kinv = (kx==0 && ky==0) ? 0.0 : -1.0/(rkx*rkx+rky*rky);
				in_out[linidx][0] *= kinv*norm_factor;
				in_out[linidx][1] *= kinv*norm_factor;
			}
		
		//this is sparta!
		in_out[0][0] = in_out[0][1] = 0;
	}
	
	void _fftw2cub(Real * out, TGrid& grid, const size_t nx, const size_t ny, const size_t ny_hat) const
	{
		vector<BlockInfo> infos = grid.getBlocksInfo();
		
		const size_t N = infos.size();
		
		const size_t mybpd[2] = {grid.getBlocksPerDimension(0), grid.getBlocksPerDimension(1)};
		
#pragma omp parallel for
		for(int i=0; i<N; ++i)
		{
			const BlockInfo info = infos[i];
			BlockType& b = *(BlockType*)infos[i].ptrBlock;
			
			const size_t offset = BlockType::sizeY*info.index[1] + (ny/2+1)*2*BlockType::sizeX*info.index[0];
			
			for(int iy=0; iy<BlockType::sizeY; iy++)
				for(int ix=0; ix<BlockType::sizeX; ix++)
				{
					const size_t src_index = TStreamer::channels*(offset + iy + 2*(ny/2+1)*ix);
					
					assert(src_index>=0 && src_index<nx*ny_hat*2); // linking error with openmp
					TStreamer::operate(&out[src_index], b(ix,iy));
				}
		}
	}
	
public:
	
	PoissonSolverScalarFFTW(const int desired_threads): FFTWBase(desired_threads), initialized(false) {	}
	
	void solve(TGrid& grid, const bool spectral=false)
	{
		const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
		const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
		
		profiler.push_start("SETUP");
		_setup(data, gsize[0], gsize[1]);
		profiler.pop_stop();
		
		profiler.push_start("CUB2FFTW");
		_cub2fftw(grid, data, gsize[0], gsize[1], gsize[1]/2+1);
		profiler.pop_stop();
		
		profiler.push_start("FFTW FORWARD");
#ifndef _SP_COMP_
		fftw_execute(fwd);
#else // _SP_COMP_
		fftwf_execute(fwd);
#endif // _SP_COMP_
		profiler.pop_stop();
		
		const Real norm_factor = 1./(gsize[0]*gsize[1]);
		const Real h = grid.getBlocksInfo().front().h_gridpoint;
		assert(1./gsize[0]==h);
		
		profiler.push_start("SOLVE");
		if(spectral)
			_solveSpectral((mycomplex *)data, gsize[0], gsize[1], gsize[1]/2+1, norm_factor, h);
		else
			_solve((mycomplex *)data, gsize[0], gsize[1], gsize[1]/2+1, norm_factor, h);
		
		profiler.pop_stop();
		
		profiler.push_start("FFTW INVERSE");
#ifndef _SP_COMP_
		fftw_execute(bwd);
#else // _SP_COMP_
		fftwf_execute(bwd);
#endif // _SP_COMP_
		profiler.pop_stop();
		
		profiler.push_start("FFTW2CUB");
		_fftw2cub(data, grid, gsize[0], gsize[1], gsize[1]/2+1);
		profiler.pop_stop();
		
		//profiler.printSummary();
	}
	
	void dispose()
	{
		if (initialized)
		{
			initialized = false;
			
#ifndef _SP_COMP_
			fftw_destroy_plan(fwd);
			fftw_destroy_plan(bwd);
			
			fftw_free(data);
#else // _SP_COMP_
			fftwf_destroy_plan(fwd);
			fftwf_destroy_plan(bwd);
			
			fftwf_free(data);
#endif // _SP_COMP_
			FFTWBase::dispose();
		}
	}
};

template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_DST : FFTWBase
{
    Profiler profiler;
    
protected:
    typedef typename TGrid::BlockType BlockType;
    
    bool initialized;
    
    myplan fwd, bwd;
    Real * data; // rhs in _setup, out in cub2fftw and fftw2cub
    
protected:
    
    virtual void _setup(Real *& rhs, const size_t nx, const size_t ny)
    {
        if (!initialized)
        {
            initialized = true;
#ifndef _SP_COMP_
            rhs = fftw_alloc_real(2*nx*(ny/2+1));
            
            fwd = fftw_plan_dft_r2c_2d(nx, ny, rhs, (mycomplex *)rhs, FFTW_MEASURE);
            bwd = fftw_plan_dft_c2r_2d(nx, ny, (mycomplex *)rhs, rhs, FFTW_MEASURE);
#else // _SP_COMP_
            rhs = fftwf_alloc_real(2*nx*(ny/2+1));
            
            fwd = fftwf_plan_dft_r2c_2d(nx, ny, rhs, (mycomplex *)rhs, FFTW_MEASURE);
            bwd = fftwf_plan_dft_c2r_2d(nx, ny, (mycomplex *)rhs, rhs, FFTW_MEASURE);
#endif // _SP_COMP_
        }
    }
    
    void _cub2fftw(TGrid& grid, Real * out, const size_t nx, const size_t ny, const size_t ny_hat) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            // indexing for fftw - [x][y] (y is the fastest running index), block part of the index
            const size_t offset = BlockType::sizeY*info.index[1] + (ny/2+1)*2*BlockType::sizeX*info.index[0];
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    // indexing for fftw - [x][y] (y is the fastest running index), local part of the index
                    const size_t dest_index = TStreamer::channels*(offset + iy + 2*(ny/2+1)*ix);
                    
                    assert(dest_index>=0 && dest_index<nx*ny_hat*2); // linking error with openmp: http://lists.apple.com/archives/xcode-users/2011/Oct/msg00252.html
                    TStreamer::operate(b(ix,iy), &out[dest_index]);
                }
        }
    }
    
    virtual void _solve(mycomplex * in_out, const size_t nx, const size_t ny, const size_t ny_hat, const Real norm_factor, const Real h)
    {
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }
        
        const Real h2 = h*h;
        const Real factor = h2*norm_factor;
        
#pragma omp parallel for
        for(int i=0; i<nx; ++i)
            for(int j = 0; j<ny_hat; ++j)
            {
                const int linidx = i*ny_hat+j;
                assert(linidx >=0 && linidx<nx*ny_hat); // linking error with openmp
                
                // based on the 5 point stencil in 1D (h^4 error)
                const Real denom = 32.*(cos(2.*M_PI*i/nx) + cos(2.*M_PI*j/ny)) - 2.*(cos(4.*M_PI*i/nx) + cos(4.*M_PI*j/ny)) - 60.;
                const Real inv_denom = (denom==0)? 0.:1./denom;
                const Real fatfactor = 12. * inv_denom * factor;
                
                // based on the 3 point stencil in 1D (h^2 error)
                //const Real denom = 2.*(cos(2.*M_PI*i/nx) + cos(2.*M_PI*j/ny)) - 4.;
                //const Real inv_denom = (denom==0)? 0.:1./denom;
                //const Real fatfactor = inv_denom * factor;
                
                // this is to check the transform only
                //const Real fatfactor = norm_factor;
                
                in_out[linidx][0] *= fatfactor;
                in_out[linidx][1] *= fatfactor;
            }
        
        //this is sparta!
        in_out[0][0] = in_out[0][1] = 0;
    }
    
    virtual void _solveSpectral(mycomplex * in_out, const size_t nx, const size_t ny, const size_t ny_hat, const Real norm_factor, const Real h)
    {
        if (TStreamer::channels != 1)
        {
            cout << "PoissonSolverScalar::PoissonSolverScalar(): Error: TStreamer::channels is " << TStreamer::channels << " and it should be 1. Aborting\n";
            abort();
        }
        
        const Real waveFactX = 2.0*M_PI/(nx*h);
        const Real waveFactY = 2.0*M_PI/(ny*h);
#pragma omp parallel for
        for(int i=0; i<nx; ++i)
            for(int j = 0; j<ny_hat; ++j)
            {
                const int linidx = i*ny_hat+j;
                assert(linidx >=0 && linidx<nx*ny_hat); // linking error with openmp
                
                // wave number
                const int kx = (i <= nx/2) ? i : -(nx-i);
                const int ky = (j <= ny/2) ? j : -(ny-j);
                const Real rkx = kx*waveFactX;
                const Real rky = ky*waveFactY;
                const Real kinv = (kx==0 && ky==0) ? 0.0 : -1.0/(rkx*rkx+rky*rky);
                in_out[linidx][0] *= kinv*norm_factor;
                in_out[linidx][1] *= kinv*norm_factor;
            }
        
        //this is sparta!
        in_out[0][0] = in_out[0][1] = 0;
    }
    
    void _fftw2cub(Real * out, TGrid& grid, const size_t nx, const size_t ny, const size_t ny_hat) const
    {
        vector<BlockInfo> infos = grid.getBlocksInfo();
        
        const size_t N = infos.size();
        
        const size_t mybpd[2] = {grid.getBlocksPerDimension(0), grid.getBlocksPerDimension(1)};
        
#pragma omp parallel for
        for(int i=0; i<N; ++i)
        {
            const BlockInfo info = infos[i];
            BlockType& b = *(BlockType*)infos[i].ptrBlock;
            
            const size_t offset = BlockType::sizeY*info.index[1] + (ny/2+1)*2*BlockType::sizeX*info.index[0];
            
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    const size_t src_index = TStreamer::channels*(offset + iy + 2*(ny/2+1)*ix);
                    
                    assert(src_index>=0 && src_index<nx*ny_hat*2); // linking error with openmp
                    TStreamer::operate(&out[src_index], b(ix,iy));
                }
        }
    }
    
public:
    
    PoissonSolverScalarFFTW_DST(const int desired_threads): FFTWBase(desired_threads), initialized(false) {	}
    
    void solve(TGrid& grid, const bool spectral=false)
    {
        const int bs[2] = {BlockType::sizeX, BlockType::sizeY};
        const size_t gsize[2] = {grid.getBlocksPerDimension(0)*bs[0], grid.getBlocksPerDimension(1)*bs[1]};
        
        profiler.push_start("SETUP");
        _setup(data, gsize[0], gsize[1]);
        profiler.pop_stop();
        
        profiler.push_start("CUB2FFTW");
        _cub2fftw(grid, data, gsize[0], gsize[1], gsize[1]/2+1);
        profiler.pop_stop();
        
        profiler.push_start("FFTW FORWARD");
#ifndef _SP_COMP_
        fftw_execute(fwd);
#else // _SP_COMP_
        fftwf_execute(fwd);
#endif // _SP_COMP_
        profiler.pop_stop();
        
        const Real norm_factor = 1./(gsize[0]*gsize[1]);
        const Real h = grid.getBlocksInfo().front().h_gridpoint;
        assert(1./gsize[0]==h);
        
        profiler.push_start("SOLVE");
        if(spectral)
            _solveSpectral((mycomplex *)data, gsize[0], gsize[1], gsize[1]/2+1, norm_factor, h);
        else
            _solve((mycomplex *)data, gsize[0], gsize[1], gsize[1]/2+1, norm_factor, h);
        
        profiler.pop_stop();
        
        profiler.push_start("FFTW INVERSE");
#ifndef _SP_COMP_
        fftw_execute(bwd);
#else // _SP_COMP_
        fftwf_execute(bwd);
#endif // _SP_COMP_
        profiler.pop_stop();
        
        profiler.push_start("FFTW2CUB");
        _fftw2cub(data, grid, gsize[0], gsize[1], gsize[1]/2+1);
        profiler.pop_stop();
        
        //profiler.printSummary();
    }
    
    void dispose()
    {
        if (initialized)
        {
            initialized = false;
            
#ifndef _SP_COMP_
            fftw_destroy_plan(fwd);
            fftw_destroy_plan(bwd);
            
            fftw_free(data);
#else // _SP_COMP_
            fftwf_destroy_plan(fwd);
            fftwf_destroy_plan(bwd);
            
            fftwf_free(data);
#endif // _SP_COMP_
            FFTWBase::dispose();
        }
    }
};
