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
    
    FluidElement() : rho(0), u(0), v(0), chi(0), p(0), pOld(0), divU(0), tmpU(0), tmpV(0), tmp(0) {}
    
    void clear()
    {
        rho = u = v = chi = p = pOld = 0;
		tmpU = tmpV = tmp = 0;
		divU = 0;
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

typedef Grid<FluidBlock, std::allocator> FluidGrid;
typedef BlockProcessing_TBB<FluidBlock> FluidBlockProcessing;
#ifndef _PERIODIC_
typedef BlockLabBottomWall<FluidBlock, std::allocator> Lab;
#else // _PERIODIC_
typedef BlockLab<FluidBlock, std::allocator> Lab;
#endif // _PERIODIC_




#endif
