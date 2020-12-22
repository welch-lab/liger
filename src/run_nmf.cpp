#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "parameters.h"
#include "entities.h"
#include "vector.h"
#include "matrix.h"
#include "factorization.h"
#include "ls.h"
#include "mask.h"
#include "nnls.h"

// NNLS implementation based on source code written by Robert Lee 
// Original code available at https://github.com/rlee32/nmf
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

DenseMatrix r_to_cpp(const NumericMatrix& A)
{
	int n = A.nrow();
	int m = A.ncol();

	DenseMatrix X = DenseMatrix(n,m);
	for(int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			X.rowmajor[i][j] = A(i,j);
		}
	}
	X.copyRowToColumn();
	return X;
}

LowerTriangularMatrix r_to_cpp_lower_triangular(const NumericMatrix& A)
{
	int n = A.nrow();
	int m = A.ncol();
	if (n != m)
	{
		//uh oh
	}
	LowerTriangularMatrix X = LowerTriangularMatrix(n);
	int ind = 0;
	for(int i = 0; i < n; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			X.rowmajor[ind] = A(i,j);
			++ind;
		}
	}
	return X;
}

void arma_to_cpp(const mat& A, DenseMatrix& X)
{
	int n = A.n_rows;
	int m = A.n_cols;

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			X.rowmajor[i][j] = A(i,j);
		}
	}
	X.copyRowToColumn();
}

void arma_to_cpp_lower_triangular(const mat& A, LowerTriangularMatrix& X)
{
	int n = A.n_rows;
	int m = A.n_cols;
	if (n != m)
	{
		//oh no
	}
	int ind = 0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			X.rowmajor[ind] = A(i,j);
			++ind;
		}
	}
}

NumericMatrix cpp_to_r(const DenseMatrix& X)
{
	int n = X.rows;
	int m = X.cols;

	NumericMatrix A = NumericMatrix(n,m);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			A(i,j) = X.rowmajor[i][j];
		}
	}
	return A;
}

void cpp_to_arma(const DenseMatrix& X, mat& A)
{
	int n = X.rows;
	int m = X.cols;
	A.set_size(n,m);
	for(int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			A(i,j) = X.rowmajor[i][j];
		}
	}
}

// [[Rcpp::export]]
arma::mat solveNNLS(const arma::mat& C, const arma::mat& B)
{
	int max_iter_nnls = 100;
	LowerTriangularMatrix CTC(C.n_cols);
	DenseMatrix CTB = DenseMatrix(B.n_cols,C.n_cols);
	arma_to_cpp_lower_triangular(C.t() * C,CTC);
	arma_to_cpp(B.t() * C,CTB);

	DenseMatrix X(C.n_cols,B.n_cols);
	NNLS_Multiple_Input nnls_input = NNLS_Multiple_Input(&CTC,X.colmajor,CTB.rowmajor,B.n_cols,max_iter_nnls);
	nnls_multiple_cpu(nnls_input);
	X.copyColumnToRow();

	mat res;
	cpp_to_arma(X,res);
	return res;
}
