#pragma once

#include <cstdlib>

#include "vector.h"
#include "entities.h"

void forwardsubstitution(LowerTriangularMatrix& G,dtype*rhs)
{//rhs is overwritten with the answer.
	for (int row = 0; row < G.dim; ++row)
	{
		int startindex = (row*(row+1))/2;
		dtype rowsum = 0;
		for(int col = 0; col < row; ++col) rowsum += G.rowmajor[startindex + col]*rhs[col];
		int diagonalindex = startindex + row;
		rhs[row] = ( rhs[row] - rowsum ) / G.rowmajor[diagonalindex];
	}
}
void backsubstitution(LowerTriangularMatrix& G,dtype*rhs)
{//column-oriented back substitutes from the transpose of G.
	//rhs is overwritten with the answer.
	for(int col = G.dim-1; col >= 0; --col)
	{
		int startindex = (col*(col+1))/2;
		int diagonalindex = startindex+col;
		rhs[col] = rhs[col] / G.rowmajor[diagonalindex];
		for(int row=0;row < col;++row) rhs[row] = rhs[row] - rhs[col]*G.rowmajor[startindex+row];
	}
}
void matmult_ata_lowertriangular_cpu(
	LowerTriangularMatrix& C, DenseMatrix& A)
{//computes C=ATA from A.
	for(int Crow=0;Crow<C.dim;++Crow)
	{
		int Cstart = (Crow*(Crow+1))/2;
		for(int Ccol=0;Ccol<=Crow;++Ccol)
		{
			int Cindex = Cstart+Ccol;
			C.rowmajor[Cindex] = vectordot(A.colmajor[Crow],A.colmajor[Ccol],A.rows);
		}
	}
}
void matmult_ata_lowertriangular_pointers_cpu(
	LowerTriangularMatrix& C, dtype** A_colmajor,int A_rows)
{//computes C=ATA from A.
	for(int Crow=0;Crow<C.dim;++Crow)
	{
		int Cstart = (Crow*(Crow+1))/2;
		for(int Ccol=0;Ccol<=Crow;++Ccol)
		{
			int Cindex = Cstart+Ccol;
			C.rowmajor[Cindex] = vectordot(A_colmajor[Crow],A_colmajor[Ccol],A_rows);
		}
	}
}
void matvecmult_transpose_cpu(DenseMatrix& A,dtype*x,dtype*b,int factor=1)
{//Multiplies A^T * x = b.
	for(int brow=0;brow<A.cols;++brow)
	{
		b[brow] = factor*vectordot(A.colmajor[brow],x,A.rows);
	}
}
// void printv(dtype*a,dtype*b,int size)
// {
// 	for(int i=0;i<size;++i)
// 	{
// 		printf("%f, %f\n",a[i],b[i]);
// 	}
// }
void matvecmult_colmajor_cpu(DenseMatrix& A,dtype*x,dtype*b,int factor=1)
{//Multiplies A * x = b, column-major on A.
	for(int row=0;row<A.rows;++row) b[row] = 0;
	for(int col=0;col<A.cols;++col)
	{
		// printv(A.colmajor[col],b,A.rows);
		vectoradd(A.colmajor[col],b,A.rows,factor*x[col]);
		// std::cin.ignore();
	}
}
void matmatmult_colmajor_cpu(DenseMatrix& A,DenseMatrix& B,DenseMatrix& C)
{//Multiplies A * x = b, column-major on A.
	if( A.rows != C.rows or B.cols != C.cols or A.cols != B.rows ) printf("Error! Dimension mismatch in matmatmult_colmajor_cpu!\n");
	for(int bcol=0;bcol<B.cols;++bcol)
	{
		for(int crow=0;crow<C.rows;++crow) C.colmajor[bcol][crow] = 0;
		matvecmult_colmajor_cpu(A,B.colmajor[bcol],C.colmajor[bcol]);
	}
}
void randInit(dtype**colmajor,int rows,int cols,int max)
{
	for(int col=0;col<cols;++col)
	{
		for(int row=0;row<rows;++row)
		{
			colmajor[col][row] = rand()%max;
		}
	}
}

dtype FrobeniusNorm(DenseMatrix* A)
{
	dtype sum = 0;
	for(int row = 0;row<A->rows;++row)
	{
		for(int col = 0;col<A->cols;++col)
		{
			sum += A->colmajor[col][row]*A->colmajor[col][row];
		}
	}
	return sqrt(sum);
}

dtype sparsity(DenseMatrix* A)
{
	int zeros = 0;
	for(int row = 0;row<A->rows;++row)
	{
		for(int col = 0;col<A->cols;++col)
		{
			if(A->colmajor[col][row] == 0) ++zeros;
		}
	}
	return (dtype)zeros/(dtype)(A->cols*A->rows);
}
void check_colmajor(dtype**standard,dtype**check,int cols,int rows,dtype tol=1e-5)
{
	int numwrong = 0;
	for(int col=0;col<cols;++col)
	{
		for(int row = 0;row < rows;++row)
		{
			if(fabs(standard[col][row]-check[col][row])/standard[col][row] > tol )
			{
				++numwrong;
			}
		}
	}
	if(numwrong > 0)
	{
		std::cout << "Error! " << numwrong*100.0/(double)(rows*cols) << " percent wrong." << std::endl;
	}
	else
	{
		std::cout << "Passed the matrix match check." << std::endl;
	}
}
void copy_colmajor(dtype**from,dtype**to,int cols,int rows)
{
	for(int col=0;col<cols;++col)
	{
		for(int row = 0;row < rows;++row)
		{
			to[col][row] = from[col][row];
		}
	}
}