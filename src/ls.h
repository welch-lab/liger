#pragma once

#include "matrix.h"
#include "factorization.h"
#include <map>
#include <string>

void normal_equations_cpu(DenseMatrix& C,dtype*x,dtype*b)
{
	//Algorithm from 5.3.2 (pg. 262) in Golub and Van Loan.
	LowerTriangularMatrix G = LowerTriangularMatrix(C.cols);
	LowerTriangularMatrix CTC = LowerTriangularMatrix(C.cols);

	matmult_ata_lowertriangular_cpu(CTC,C);
	matvecmult_transpose_cpu(C,b,x);
	cholesky_lowertriangular_cpu(G,CTC);
	
	forwardsubstitution(G,x);
	backsubstitution(G,x);

}

void normal_equations_precomputedCholesky_cpu(LowerTriangularMatrix& G,dtype*x)
{//x needs to be initialized with appropriately masked CTb.
	forwardsubstitution(G,x);
	backsubstitution(G,x);
}

void normal_equations_precomputedCholesky_cpu_test(LowerTriangularMatrix& G,dtype*x)
{//x needs to be initialized with appropriately masked CTb.
	forwardsubstitution(G,x);
	backsubstitution(G,x);
}


