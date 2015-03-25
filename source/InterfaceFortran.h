//
//  InterfaceFortran.h
//  CubismUP_2D
//
//  Created by Christian Conti on 2/6/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_InterfaceFortran_h
#define CubismUP_2D_InterfaceFortran_h

void cofx(int ix, float& ax, float& bx, float& cx)
{
	ax = 1;
	bx = 0;
	cx = 0;
}

// fishpack interface
extern "C"
{
	// remember:
	// transposition (column vs row major)
	// indices (0 vs 1 start index)
	void sepx4_(int * IORDER, 	// order of approximation (2 or 4)
				float * A,		// start of x interval
				float * B,		// end of x interval
				int * M,		// M+1 gridpoints in x
				int * MBDCND,	// BC - 0: periodic in x
				//		1: dirichlet in x
				//		2: dirichlet in a, mixed in b
				//		3: mixed bc in x
				//		4: dirichlet in b, mixed in a
				float * BDA,		// 1D array of size N+1 containing the BC on A for choices 3,4
				float * ALPHA,	// BC stuff
				float * BDB,		// 1D array of size N+1 containing the BC on B for choices 2,3
				float * BETA,	// BC stuff
				float * C,		// start of y interval
				float * D,		// end of y interval
				int * N,			// N+1 gridpoints in y
				int * NBDCND,	// BC in y direction
				float * BDC,		// see BDA
				float * BDD,		// see BDB
				void (*COFX)(int, float&, float&, float&),	// some subprogram - function pointer (?)
				float * GRHS,	// RHS 2D array
				float * USOL,	// 2D boundary solution array (in/out)
				int * IDMN,		// row dimension of GRHS,USOL
				float * W,		// 1D working space (in/out)
				float * PERTRB,	// const value subtracted to ensure existance of solution (out)
				int * IERROR);	// integer error flag (out)
}

// fishpack interface v2
extern "C"
{
	void hwscrt_(float * A,			// start of x interval
				 float * B,			// end of x interval
				 int * M,			// M+1 gridpoints in x
				 int * MBDCND,		// BC - 0: periodic in x
				 					//		1: dirichlet in x
				 					//		2: dirichlet in a, mixed in b
				 					//		3: mixed bc in x
				 					//		4: dirichlet in b, mixed in a
				 float * BDA,		// 1D array of size N+1 containing the BC on A for choices 3,4
				 float * BDB,		// 1D array of size N+1 containing the BC on B for choices 2,3
				 float * C,			// start of y interval
				 float * D,			// end of y interval
				 int * N,			// N+1 gridpoints in y
				 int * NBDCND,		// BC in y direction
				 float * BDC,		// see BDA
				 float * BDD,		// see BDB
				 float * ELMBDA, 	// constant in the helmhotz equation multiplying the unknown function
				 float * F,			// RHS 2D array, solution on output
				 int * IDIMF,		// row dimension of GRHS,USOL
				 float * PERTRB,	// const value subtracted to ensure existance of solution (out)
				 int * IERROR);		// integer error flag (out)
}

#endif
