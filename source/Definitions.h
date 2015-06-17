//
//  DataStructures.h
//  CubismUP_2D
//
//  Created by Christian Conti on 1/7/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_DataStructures_h
#define CubismUP_2D_DataStructures_h

#include "common.h"
#include "Layer.h"
#include "LayerToVTK.h"
#include "BoundaryConditions.h"

#ifndef _BS_
#define _BS_ 32
#endif // _BS_

struct FluidElement
{
    Real rho, u, v, chi, p, pOld;
	Real tmpU, tmpV, tmp;
	Real divU;
	Real x, y;
    
    FluidElement() : rho(0), u(0), v(0), chi(0), p(0), pOld(0), divU(0), tmpU(0), tmpV(0), tmp(0), x(0), y(0) {}
    
    void clear()
    {
        rho = u = v = chi = p = pOld = 0;
		tmpU = tmpV = tmp = 0;
		divU = 0;
		x = y = 0;
    }
};


struct FluidVTKStreamer
{
	static const int channels = 8;
	
	void operate(FluidElement input, Real output[8])
	{
		output[0] = input.rho;
		output[1] = input.u;
		output[2] = input.v;
		output[3] = input.p;
		output[4] = input.chi;
		output[5] = input.divU;
		output[6] = input.pOld;
		output[7] = input.tmp;
	}
};

// this is used for serialization - important that ALL the quantities are streamed
struct StreamerGridPoint
{
	static const int channels = 10;
	
	void operate(const FluidElement& input, Real output[10]) const
	{
		abort();
		output[0] = input.rho;
		output[1] = input.u;
		output[2] = input.v;
		output[3] = input.chi;
		output[4] = input.p;
		output[5] = input.pOld;
		output[6] = input.tmpU;
		output[7] = input.tmpV;
		output[8] = input.tmp;
		output[9] = input.divU;
	}
	
	void operate(const Real input[10], FluidElement& output) const
	{
		abort();
		output.rho  = input[0];
		output.u    = input[1];
		output.v    = input[2];
		output.chi  = input[3];
		output.p    = input[4];
		output.pOld = input[5];
		output.tmpU = input[6];
		output.tmpV = input[7];
		output.tmp  = input[8];
		output.divU = input[9];
	}
};

struct StreamerGridPointASCII
{
	void operate(const FluidElement& input, ofstream& output) const
	{
		output << input.rho << " " << input.u << " " << input.v << " " << input.chi << " " << input.p << " " << input.pOld << " " << input.tmpU << " " << input.tmpV << " " << input.tmp << " " << input.divU;
	}
	
	void operate(ifstream& input, FluidElement& output) const
	{
		input >> output.rho;
		input >> output.u;
		input >> output.v;
		input >> output.chi;
		input >> output.p;
		input >> output.pOld;
		input >> output.tmpU;
		input >> output.tmpV;
		input >> output.tmp;
		input >> output.divU;
	}
};

struct StreamerDiv
{
	static const int channels = 1;
	static void operate(const FluidElement& input, Real output[1])
	{
		output[0] = input.divU;
	}
	
	static void operate(const Real input[1], FluidElement& output)
	{
		output.divU = input[0];
	}
};

struct FluidBlock
{
    //these identifiers are required by cubism!
    static const int sizeX = _BS_;
    static const int sizeY = _BS_;
    static const int sizeZ = 1;
    typedef FluidElement ElementType;
    FluidElement data[1][sizeY][sizeX];
    
    //required from Grid.h
    void clear()
    {
        FluidElement * entry = &data[0][0][0];
        const int N = sizeX*sizeY;
        
        for(int i=0; i<N; ++i)
            entry[i].clear();
    }
    
    FluidElement& operator()(int ix, int iy=0, int iz=0)
    {
        assert(ix>=0); assert(ix<sizeX);
        assert(iy>=0); assert(iy<sizeY);
        
        return data[0][iy][ix];
	}
	
	template <typename Streamer>
	inline void Write(ofstream& output, Streamer streamer) const
	{
		for(int iy=0; iy<sizeY; iy++)
			for(int ix=0; ix<sizeX; ix++)
				streamer.operate(data[0][iy][ix], output);
	}
	
	template <typename Streamer>
	inline void Read(ifstream& input, Streamer streamer)
	{
		for(int iy=0; iy<sizeY; iy++)
			for(int ix=0; ix<sizeX; ix++)
				streamer.operate(input, data[0][iy][ix]);
	}
};

template <> inline void FluidBlock::Write<StreamerGridPoint>(ofstream& output, StreamerGridPoint streamer) const
{
	output.write((const char *)&data[0][0][0], sizeof(FluidElement)*sizeX*sizeY);
}

template <> inline void FluidBlock::Read<StreamerGridPoint>(ifstream& input, StreamerGridPoint streamer)
{
	input.read((char *)&data[0][0][0], sizeof(FluidElement)*sizeX*sizeY);
}



struct StreamerSerialization
{
	static const int NCHANNELS = 10;
	
	FluidBlock& ref;
	
	StreamerSerialization(FluidBlock& b): ref(b) {}
	
	void operate(const int ix, const int iy, const int iz, Real output[10]) const
	{
		const FluidElement& input = ref.data[iz][iy][ix];
		
		output[0] = input.rho;
		output[1] = input.u;
		output[2] = input.v;
		output[3] = input.chi;
		output[4] = input.p;
		output[5] = input.pOld;
		output[6] = input.tmpU;
		output[7] = input.tmpV;
		output[8] = input.tmp;
		output[9] = input.divU;
	}
	
	void operate(const Real input[10], const int ix, const int iy, const int iz) const
	{
		FluidElement& output = ref.data[iz][iy][ix];
		
		output.rho  = input[0];
		output.u    = input[1];
		output.v    = input[2];
		output.chi  = input[3];
		output.p    = input[4];
		output.pOld = input[5];
		output.tmpU = input[6];
		output.tmpV = input[7];
		output.tmp  = input[8];
		output.divU = input[9];
	}
	
	void operate(const int ix, const int iy, const int iz, Real *ovalue, const int field) const
	{
		const FluidElement& input = ref.data[iz][iy][ix];
		
		switch(field) {
			case 0: *ovalue = input.rho; break;
			case 1: *ovalue = input.u; break;
			case 2: *ovalue = input.v; break;
			case 3: *ovalue = input.chi; break;
			case 4: *ovalue = input.p; break;
			case 5: *ovalue = input.pOld; break;
			case 6: *ovalue = input.tmpU; break;
			case 7: *ovalue = input.tmpV; break;
			case 8: *ovalue = input.tmp; break;
			case 9: *ovalue = input.divU; break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}
	
	void operate(const Real ivalue, const int ix, const int iy, const int iz, const int field) const
	{
		FluidElement& output = ref.data[iz][iy][ix];
		
		switch(field) {
			case 0:  output.rho  = ivalue; break;
			case 1:  output.u    = ivalue; break;
			case 2:  output.v    = ivalue; break;
			case 3:  output.chi  = ivalue; break;
			case 4:  output.p    = ivalue; break;
			case 5:  output.pOld = ivalue; break;
			case 6:  output.tmpU = ivalue; break;
			case 7:  output.tmpV = ivalue; break;
			case 8:  output.tmp  = ivalue; break;
			case 9:  output.divU = ivalue; break;
			default: throw std::invalid_argument("unknown field!"); break;
		}
	}
	
	static const char * getAttributeName() { return "Tensor"; }
};


template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabDirichlet : public BlockLab<BlockType,allocator>
{
	typedef typename BlockType::ElementType ElementTypeBlock;
	
public:
	BlockLabDirichlet(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{
		BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
		
		ElementTypeBlock p;
		
		p.rho = 1;
		p.u   = 0;
		p.v   = 0;
		p.tmp = 0;
		if (info.index[0]==0)          bc.template applyBC_dirichlet<0,0>(p);
		if (info.index[0]==this->NX-1) bc.template applyBC_dirichlet<0,1>(p);
		if (info.index[1]==0)		   bc.template applyBC_dirichlet<1,0>(p);
		if (info.index[1]==this->NY-1) bc.template applyBC_dirichlet<1,1>(p);
	}
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabNeumann : public BlockLab<BlockType,allocator>
{
	typedef typename BlockType::ElementType ElementTypeBlock;
	
public:
	BlockLabNeumann(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{
		BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
		
		ElementTypeBlock p;
		
		if (info.index[0]==0)          bc.template applyBC_neumann<0,0>();
		if (info.index[0]==this->NX-1) bc.template applyBC_neumann<0,1>();
		if (info.index[1]==0)		   bc.template applyBC_neumann<1,0>();
		if (info.index[1]==this->NY-1) bc.template applyBC_neumann<1,1>();
	}
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabBottomWall : public BlockLab<BlockType,allocator>
{
	typedef typename BlockType::ElementType ElementTypeBlock;
	
public:
	BlockLabBottomWall(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{
		BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
		
		// keep periodicity in x direction
		if (info.index[1]==0)		   bc.template applyBC_mixedBottom<1,0>();
		if (info.index[1]==this->NY-1) bc.template applyBC_mixedTop<1,1>();
		//if (info.index[1]==this->NY-1) bc.template applyBC_mixedBottom<1,1>();
	}
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabPipe : public BlockLab<BlockType,allocator>
{
	typedef typename BlockType::ElementType ElementTypeBlock;
	
public:
	BlockLabPipe(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{
		BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
		
		if (info.index[1]==0)		   bc.template applyBC_mixedBottom<1,0>();
		if (info.index[1]==this->NY-1) bc.template applyBC_mixedBottom<1,1>();
	}
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabBox : public BlockLab<BlockType,allocator>
{
	typedef typename BlockType::ElementType ElementTypeBlock;
	
public:
	BlockLabBox(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{
		BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
		
		if (info.index[0]==0)		   bc.template applyBC_mixedBottom<0,0>();
		if (info.index[0]==this->NX-1) bc.template applyBC_mixedBottom<0,1>();
		if (info.index[1]==0)		   bc.template applyBC_mixedBottom<1,0>();
		if (info.index[1]==this->NY-1) bc.template applyBC_mixedTop<1,1>();
	}
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabVortex : public BlockLab<BlockType,allocator>
{
	typedef typename BlockType::ElementType ElementTypeBlock;
	
public:
	BlockLabVortex(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{
		BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
		
		if (info.index[0]==0)		   bc.template applyBC_vortex<0,0>(info);
		if (info.index[0]==this->NX-1) bc.template applyBC_vortex<0,1>(info);
		if (info.index[1]==0)		   bc.template applyBC_vortex<1,0>(info);
		if (info.index[1]==this->NY-1) bc.template applyBC_vortex<1,1>(info);
	}
};

typedef Grid<FluidBlock, std::allocator> FluidGrid;

#ifdef _MIXED_
typedef BlockLabBottomWall<FluidBlock, std::allocator> Lab;
#endif // _MIXED_

#ifdef _PERIODIC_
typedef BlockLab<FluidBlock, std::allocator> Lab;
#endif // _PERIODIC_

#ifdef _VORTEX_
typedef BlockLabVortex<FluidBlock, std::allocator> Lab;
#endif // _VORTEX_

#ifdef _PIPE_
typedef BlockLabPipe<FluidBlock, std::allocator> Lab;
#endif // _PIPE_




#endif
