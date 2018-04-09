#include "debug.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;
using std::cout;
using std::endl;

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
		cout << "Warning: Asymmetric matrix in r_to_cpp_lower_triangular" << endl;
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

	for(int i = 0; i < n; ++i)
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
		cout << "Warning: Non-square matrix in r_to_cpp_lower_triangular" << endl;
	}
	int ind = 0;
	for(int i = 0; i < n; ++i)
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
	for(int i = 0; i < n; ++i)
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

double norm_sq(const DenseMatrix& X)
{
	std::cout << "norm_sq" << std::endl;
	int n = X.rows;
	int m = X.cols;
	double norm_val = 0;
	for(int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			norm_val += X.rowmajor[i][j]*X.rowmajor[i][j];
		}
	}
	return norm_val;
}

double calculate_objective(const DenseMatrix& E1, const DenseMatrix& H1, const DenseMatrix& W, const DenseMatrix& V1, const DenseMatrix& E2, const DenseMatrix& H2, const DenseMatrix& V2, double lambda)
{
	return(norm_sq(E1-H1*(W+V1))+norm_sq(E2-H2*(W+V2))+lambda*norm_sq(H1*V1)+lambda*norm_sq(H2*V2));
}

double calculate_objective(const DenseMatrix& A, const DenseMatrix& W, const DenseMatrix& H)
{
	mat A_arma, W_arma, H_arma;
	cpp_to_arma(A,A_arma);
	cpp_to_arma(W,W_arma);
	cpp_to_arma(H,H_arma);
	double frob_norm = norm(A_arma - (W_arma * H_arma),"fro");
	return (frob_norm*frob_norm);
}

// [[Rcpp::export]]
arma::mat solve_nnls(const arma::mat& C, const arma::mat& B)
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

List nmf_als(NumericMatrix A, mat W_in, mat H_in, int k, double thresh = 0.0001)
{
	int n = A.nrow();	
	int m = A.ncol();
	DenseMatrix X = r_to_cpp(A);
	DenseMatrix W = DenseMatrix(n,k);
	arma_to_cpp(W_in,W);
	DenseMatrix H = DenseMatrix(k,m);
	arma_to_cpp(H_in,H);
	NMF_Input input = NMF_Input(&W,&H,&X,1,100);
	NMF_State state = NMF_State(input.m,input.k,input.n);
	NNLS_Multiple_Input nnls_input_1 = NNLS_Multiple_Input(state.HHT,input.W->rowmajor,state.HAT,input.m,input.max_iter_nnls);
	NNLS_Multiple_Input nnls_input_2 = NNLS_Multiple_Input(state.WTW,input.H->colmajor,state.WTA,input.n,input.max_iter_nnls);
	
	//initializeMatrix(input.H);
	int iterations = 0;
	int iterations_nnls = 0;
	std::cout << "Obj0" << std::endl;
	double obj = calculate_objective(X,W,H);
	std::cout << "After obj0" << std::endl;
	double obj0;
	double delta = 1.0;
	while(delta > thresh && iterations < input.max_iter_nmf)//TODO: implement solution-sensitive stopping criterion
	{
		obj0 = obj;
		std::cout << obj << std::endl;
		input.H->copyColumnToRow();
		
		cout << "W:\n" << *input.W << endl;
		cout << "H:\n" << *input.H << endl;
		
		matmult_ata_lowertriangular_pointers_cpu(*state.HHT,input.H->rowmajor,input.H->cols);
		for(int i=0;i<input.m;++i) matvecmult_colmajor_cpu(*input.H,input.A->rowmajor[i],state.HAT[i]);
		iterations_nnls += nnls_multiple_cpu(nnls_input_1);
		input.W->copyRowToColumn();
		
		cout << "HHT:\n" << *state.HHT << endl;
		cout << "HAT:\n";
		for(int i=0;i<n;++i)
		{
			for(int j=0;j<k;++j)
			{
			 	cout << state.HAT[i][j] << " ";
			}
			cout << endl;
		}
		cout << "W:\n" << *input.W << endl;
		
		matmult_ata_lowertriangular_cpu(*state.WTW,*input.W);
		for(int i=0;i<input.n;++i) matvecmult_transpose_cpu(*input.W,input.A->colmajor[i],state.WTA[i]);
		iterations_nnls += nnls_multiple_cpu(nnls_input_2);
		
		cout << "WTW:\n" << *state.WTW << endl;
		cout << "WTA:\n";
		for(int i=0;i<m;++i)
		{
			for(int j=0;j<k;++j)
			{
			 	cout << state.WTA[i][j] << " ";
			}
			cout << endl;
		}
		cout << "H:\n" << *input.H << endl;
		
		List nmf_res;
		return nmf_res;
		
		++iterations;
		obj = calculate_objective(X,W,H);
		delta = (obj0-obj)/(obj0);
	}
		
	NumericMatrix W_ret = cpp_to_r(W);
	NumericMatrix H_ret = cpp_to_r(X);
	List nmf_res;
	nmf_res["W"]=W_ret;
	nmf_res["H"]=H_ret;
	return nmf_res;
}

/*List iNMF_als(NumericMatrix E1, NumericMatrix E2, int k, double lambda=5.0, double thresh=0.0001, int max_iters, nrep=20)
{
	n1 = E1.nrow()
	n2 = E2.nrow()
	g = E1.ncol()
	
	NMF_State state = NMF_State(input.m,input.k,input.n);
	NNLS_Multiple_Input nnls_input_H1 = NNLS_Multiple_Input(state.HHT,input.W->rowmajor,state.HAT,input.m,input.max_iter_nnls);
	NNLS_Multiple_Input nnls_input_H22 = NNLS_Multiple_Input(state.WTW,input.H->colmajor,state.WTA,input.n,input.max_iter_nnls);

	// initializeMatrix(input.H,averageMatrix(input.A));
	initializeMatrix(input.H);
	int iterations = 0;
	int iterations_nnls = 0;
	obj0 = calculate_objective(E1,E2,H1,H2,W,V1,V2);
	while(delta < thresh)//TODO: implement solution-sensitive stopping criterion
	{
		input.H->copyColumnToRow();
		matmult_ata_lowertriangular_pointers_cpu(*state.HHT,input.H->rowmajor,input.H->cols);
		for(int i=0;i<input.m;++i) matvecmult_colmajor_cpu(*input.H,input.A->rowmajor[i],state.HAT[i]);
		iterations_nnls += nnls_multiple_cpu(nnls_input_1);

		input.W->copyRowToColumn();
		matmult_ata_lowertriangular_cpu(*state.WTW,*input.W);
		for(int i=0;i<input.n;++i) matvecmult_transpose_cpu(*input.W,input.A->colmajor[i],state.WTA[i]);
		iterations_nnls += nnls_multiple_cpu(nnls_input_2);
		
		++iterations;
		obj = calculate_objective(E1,E2,H1,H2,W,V1,V2);
		delta = abs(obj-obj0)/((obj+obj0)/2);
	}
	
	W1_m = matrix(abs(rnorm(g * k)), k, g)
	V1_m = matrix(abs(rnorm(g * k)), k, g)
	V2_m = matrix(abs(rnorm(g * k)), k, g)
	H1_m = matrix(abs(rnorm(n1 * k)), n1, k)
	H2_m = matrix(abs(rnorm(n2 * k)), n2, k)
	best_obj = Inf
	run_stats = matrix(0,nrow=nrep,ncol=2)
	for (i in 1:nrep)
	{
	start_time <- Sys.time()
	W1 = matrix(abs(runif(g * k,0,2)), k, g)
	V1 = matrix(abs(runif(g * k,0,2)), k, g)
	V2 = matrix(abs(runif(g * k,0,2)), k, g)
	H1 = matrix(abs(runif(n1 * k,0,2)), n1, k)
	H2 = matrix(abs(runif(n2 * k,0,2)), n2, k)
	delta = 1
	iters = 0
	#pb = txtProgressBar(min=0,max=max_iters,style=3)
	sqrt_lambda = sqrt(lambda)
	obj0 = norm(E1-H1%*%(W1+V1),"F")^2+norm(E2-H2%*%(W1+V2),"F")^2+lambda*norm(H1%*%V1,"F")^2+lambda*norm(H2%*%V2,"F")^2
	start_obj = obj0
	while(delta > thresh & iters < max_iters)
	{
	  H1 = t(fcnnls(x=rbind(t(W1)+t(V1),sqrt_lambda*t(V1)),y=rbind(t(E1),matrix(0,nrow=g,ncol=n1)))$x)
	  H2 = t(fcnnls(x=rbind(t(W1)+t(V2),sqrt_lambda*t(V2)),y=rbind(t(E2),matrix(0,nrow=g,ncol=n2)))$x)
	  V1 = fcnnls(x=rbind(H1,sqrt_lambda*H1),y=rbind(E1-H1%*%W1,matrix(0,nrow=n1,ncol=g)))$x
	  V2 = fcnnls(x=rbind(H2,sqrt_lambda*H2),y=rbind(E2-H2%*%W1,matrix(0,nrow=n2,ncol=g)))$x
	  W1 = fcnnls(x=rbind(H1,H2),y=rbind(E1-H1%*%V1,E2-H2%*%V2),verbose=T)$x
	  obj = norm(E1-H1%*%(W1+V1),"F")^2+norm(E2-H2%*%(W1+V2),"F")^2+lambda*norm(H1%*%V1,"F")^2+lambda*norm(H2%*%V2,"F")^2
	  print(obj)
	  delta = abs(obj0-obj)/(mean(obj0,obj))
	  #delta = abs(obj - obj0)/(start_obj-obj)
	  obj0 = obj
	  iters = iters + 1
	  #setTxtProgressBar(pb,iters)
	}
	#setTxtProgressBar(pb,max_iters)
	if (obj<best_obj)
	{
	  W1_m = W1
	  H1_m = H1
	  H2_m = H2
	  V1_m = V1
	  V2_m = V2
	  best_obj = obj
	}
	end_time <- Sys.time()
	run_stats[i,1]=as.double(difftime(end_time,start_time,units="secs"))
	run_stats[i,2]=iters
	cat("\nConverged in",run_stats[i,1],"seconds,",iters,"iterations. Objective:",obj,"\n")
	}
	cat("\nBest objective:",best_obj,"\n")
	return(list(H1_m,H2_m,W1_m,V1_m,V2_m,run_stats))
}*/