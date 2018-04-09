#pragma once

#include "matrix.h"
#include "factorization.h"
#include <map>
#include <string>

void printv(dtype*a,int size)
{
	for(int i=0;i<size;++i)
	{
		printf("%f,",a[i]);
	}
	printf("\n");
}
void normal_equations_cpu(DenseMatrix& C,dtype*x,dtype*b)
{
	//Algorithm from 5.3.2 (pg. 262) in Golub and Van Loan.
	LowerTriangularMatrix G = LowerTriangularMatrix(C.cols);
	LowerTriangularMatrix CTC = LowerTriangularMatrix(C.cols);

	matmult_ata_lowertriangular_cpu(CTC,C);
	matvecmult_transpose_cpu(C,b,x);
	cholesky_lowertriangular_cpu(G,CTC);
	// printv(G.rowmajor,G.totalsize);
	// printf("start\n");
	// printv(x,G.dim);
	forwardsubstitution(G,x);
	// printv(x,G.dim);
	backsubstitution(G,x);
	// printv(x,G.dim);

}

void normal_equations_precomputedCholesky_cpu(LowerTriangularMatrix& G,dtype*x)
{//x needs to be initialized with appropriately masked CTb.
	// printf("start\n");
	// printv(x,G.dim);
	forwardsubstitution(G,x);
	// printv(x,G.dim);
	backsubstitution(G,x);
	// printv(x,G.dim);
}

void normal_equations_precomputedCholesky_cpu_test(LowerTriangularMatrix& G,dtype*x)
{//x needs to be initialized with appropriately masked CTb.
	// printf("start\n");
	// printv(x,G.dim);
	forwardsubstitution(G,x);
	// printv(x,G.dim);
	backsubstitution(G,x);
	// printv(x,G.dim);
}


