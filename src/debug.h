#pragma once

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
#include "nmf.h"
#include "reader.h"

void checkvector(dtype*a,dtype*b,int size,std::string testname,dtype tol=0)
{
	bool match = true;
	for(int i=0;i<size;++i)
	{
		if ( fabs(a[i]-b[i]) > tol )
		{
			match = false;
			break;
		}
	}
	std::string printstring(testname);
	if(match)
		printstring = "Passed "+printstring;
	else
		printstring = "Failed "+printstring;
	std::cout << printstring << "." << std::endl;
}
void checkMatrix(dtype** M,dtype** M_expected,int rows,int cols,std::string testname,dtype tol=0)
{
	bool match = true;
	for(int row=0;row<rows;++row)
	{
		for(int col=0;col<cols;++col)
		{
			if ( fabs(M[row][col]-M_expected[row][col]) > tol )
			{
				match = false;
				break;
			}
		}
		if( not match ) break;
	}
	std::string printstring(testname);
	if(match)
		printstring = "Passed "+printstring;
	else
		printstring = "Failed "+printstring;
	std::cout << printstring << "." << std::endl;
}
void printrowmajor(dtype** M,int rows,int cols)
{
	for(int row=0;row<rows;++row)
	{
		for(int col=0;col<cols;++col)
		{
			printf("%f,",M[row][col]);
		}
		printf("\n");
	}
}
void printcolmajor(dtype** M,int rows,int cols)
{
	for(int row=0;row<rows;++row)
	{
		for(int col=0;col<cols;++col)
		{
			printf("%f,",M[col][row]);
		}
		printf("\n");
	}
}
dtype Fnorm(DenseMatrix& W,DenseMatrix& H,DenseMatrix& A)
{
	if (A.rows != W.rows or A.cols != H.cols or W.cols != H.rows)
	{
		printf("Warning! Input matrix dimensions do not match in the Frobenius norm calculation.\n");
	}
	dtype dsum = 0;
	for(int col=0;col<A.cols;++col)
	{	
		for(int row=0;row<A.rows;++row)
		{
			dtype aelem = A.rowmajor[row][col];
			//dtype sum = 0;
			dtype* a = W.rowmajor[row];
			dtype* b = H.colmajor[col];
			dtype sum = vectordot(a,b,W.cols);
			// printf("Computed W*H (%d,%d): %f\n",row,col,sum);
			//for(int i=0;i<W.cols;++i) {printf("(%d,%f,%f)\t",i,a[i],b[i]);sum += a[i]*b[i];}
			// printf("\n");
			dtype diff = sum-aelem;
			// printf("col,row,diff,aelem: %d,%d,%f,%f\n",col,row,diff,aelem);
			dsum += diff*diff;
		}
	}
	return sqrt(dsum);
}
void printvectors(dtype*a,dtype*b,int size)
{
	for(int i=0;i<size;++i)
	{
		printf("%f, %f\n",a[i],b[i]);
	}
}


void test_substitution_cpu()
{
	LowerTriangularMatrix G = LowerTriangularMatrix(4);

	G.rowmajor[0] = 3;
	G.rowmajor[1] = -1;
	G.rowmajor[2] = 1;
	G.rowmajor[3] = 3;
	G.rowmajor[4] = -2;
	G.rowmajor[5] = -1;
	G.rowmajor[6] = 1;
	G.rowmajor[7] = -2;
	G.rowmajor[8] = 6;
	G.rowmajor[9] = 2;
	dtype rhs[4] = {5, 6, 4, 2};
	dtype ans[4] = {5./3., 23./3., -43./3., 305./6.};
	forwardsubstitution(G,rhs);
	checkvector(rhs,ans,4,"forwardsubstitution");

	G.rowmajor[0] = 4;
	G.rowmajor[1] = -1;
	G.rowmajor[2] = -2;
	G.rowmajor[3] = 2;
	G.rowmajor[4] = 7;
	G.rowmajor[5] = 6;
	G.rowmajor[6] = 3;
	G.rowmajor[7] = -4;
	G.rowmajor[8] = 5;
	G.rowmajor[9] = 3;
	rhs[0] = 20;rhs[1] = -7;rhs[2] = 4;rhs[3] = 6;
	ans[0] = 3;ans[1] = -4;ans[2] = -1;ans[3] = 2;
	backsubstitution(G,rhs);
	checkvector(rhs,ans,4,"backsubstitution");
}

void test_matmult_ata_lowertriangular_cpu()
{
	DenseMatrix A = DenseMatrix(4,4);
	A.colmajor[0][0] = 2;A.colmajor[1][0] = 7;A.colmajor[2][0] = 9;A.colmajor[3][0] = 7;
	A.colmajor[0][1] = 8;A.colmajor[1][1] = 6;A.colmajor[2][1] = 1;A.colmajor[3][1] = 2;
	A.colmajor[0][2] = 0;A.colmajor[1][2] = 8;A.colmajor[2][2] = 6;A.colmajor[3][2] = 5;
	A.colmajor[0][3] = 3;A.colmajor[1][3] = 7;A.colmajor[2][3] = 7;A.colmajor[3][3] = 2;
	LowerTriangularMatrix C = LowerTriangularMatrix(4);
	LowerTriangularMatrix C_expected = LowerTriangularMatrix(4);
	matmult_ata_lowertriangular_cpu(C,A);
	C_expected.rowmajor[0] = 77;
	C_expected.rowmajor[1] = 83;C_expected.rowmajor[2] = 198;
	C_expected.rowmajor[3] = 47;C_expected.rowmajor[4] = 166;C_expected.rowmajor[5] = 167;
	C_expected.rowmajor[6] = 36;C_expected.rowmajor[7] = 115;C_expected.rowmajor[8] = 109;C_expected.rowmajor[9] = 82;
	checkvector(C.rowmajor,C_expected.rowmajor,C_expected.totalsize,"test_matmult_ata_lowertriangular_cpu");
}

void test_matvecmult_transpose_cpu()
{
	DenseMatrix A = DenseMatrix(4,4);
	A.colmajor[0][0] = 9;A.colmajor[1][0] = 3;A.colmajor[2][0] = 6;
	A.colmajor[0][1] = 2;A.colmajor[1][1] = 6;A.colmajor[2][1] = 5;
	A.colmajor[0][2] = 3;A.colmajor[1][2] = 5;A.colmajor[2][2] = 4;
	A.colmajor[0][3] = 4;A.colmajor[1][3] = 3;A.colmajor[2][3] = 3;
	dtype x[4];
	x[0] = 2;
	x[1] = 9;
	x[2] = 7;
	x[3] = 3;
	dtype b_expected[3];
	b_expected[0] = 69;
	b_expected[1] = 104;
	b_expected[2] = 94;
	dtype b[3];
	matvecmult_transpose_cpu(A,x,b);
	checkvector(b,b_expected,3,"test_matvecmult_transpose_cpu");
}


void test_cholesky_lowertriangular_cpu()
{
	LowerTriangularMatrix C = LowerTriangularMatrix(4);
	LowerTriangularMatrix G = LowerTriangularMatrix(4);
	LowerTriangularMatrix G_expected = LowerTriangularMatrix(4);

	C.rowmajor[0] = 61;
	C.rowmajor[1] = 5;
	C.rowmajor[2] = 31;
	C.rowmajor[3] = 10;
	C.rowmajor[4] = 5;
	C.rowmajor[5] = 21;
	C.rowmajor[6] = 1;
	C.rowmajor[7] = 4;
	C.rowmajor[8] = 2;
	C.rowmajor[9] = 81;

	G_expected.rowmajor[0] = 7.8102497;
	G_expected.rowmajor[1] = 0.6401844;
	G_expected.rowmajor[2] = 5.5308375;
	G_expected.rowmajor[3] = 1.2803688;
	G_expected.rowmajor[4] = 0.7558219;
	G_expected.rowmajor[5] = 4.3346729;
	G_expected.rowmajor[6] = 0.1280369;
	G_expected.rowmajor[7] = 0.7083977;
	G_expected.rowmajor[8] = 0.3000556;
	G_expected.rowmajor[9] = 8.9661444;

	cholesky_lowertriangular_cpu(G,C);
	checkvector(G.rowmajor,G_expected.rowmajor,
		G_expected.totalsize,"cholesky_lowertriangular_cpu",1e-5);

}


void test_normal_equations_cpu()
{
	DenseMatrix A = DenseMatrix(4,4);
	A.colmajor[0][0] = 2;A.colmajor[1][0] = 7;A.colmajor[2][0] = 9;A.colmajor[3][0] = 7;
	A.colmajor[0][1] = 8;A.colmajor[1][1] = 6;A.colmajor[2][1] = 1;A.colmajor[3][1] = 2;
	A.colmajor[0][2] = 0;A.colmajor[1][2] = 8;A.colmajor[2][2] = 6;A.colmajor[3][2] = 5;
	A.colmajor[0][3] = 3;A.colmajor[1][3] = 7;A.colmajor[2][3] = 7;A.colmajor[3][3] = 2;
	dtype b[4];
	b[0] = 2;
	b[1] = 2;
	b[2] = 9;
	b[3] = 7;
	dtype x[4];
	dtype x_expected[4];
	x_expected[0] = -0.9612263;
	x_expected[1] = 2.0225428;
	x_expected[2] = -0.3047791;
	x_expected[3] = -1.0703336;

	normal_equations_cpu(A,x,b);

	checkvector(x,x_expected,
		4,"normal_equations_cpu",1e-5);

}


void test_markInfeasible()
{
	Mask xmask = Mask(8,false);
	xmask.value[0] = true;
	xmask.value[3] = true;
	xmask.value[5] = true;
	xmask.value[6] = true;
	Mask infeasiblemask = Mask(8,false);
	dtype x_masked[4] = {-1, 2, 3, -6};
	dtype y_masked[4] = {-1, -2, 3, 4};

	int infeasible = markInfeasible(infeasiblemask,x_masked,y_masked,xmask);
	if(infeasible == 4) printf("Passed test_markInfeasible.\n");
	else printf("Failed test_markInfeasible.\n");

}


void test_nnls_single_cpu()
{
	DenseMatrix A = DenseMatrix(4,4);
	A.colmajor[0][0] = 8;A.colmajor[1][0] = 6;A.colmajor[2][0] = 10;A.colmajor[3][0] = 10;
	A.colmajor[0][1] = 9;A.colmajor[1][1] = 1;A.colmajor[2][1] = 10;A.colmajor[3][1] = 5;
	A.colmajor[0][2] = 1;A.colmajor[1][2] = 3;A.colmajor[2][2] = 2;A.colmajor[3][2] = 8;
	A.colmajor[0][3] = 9;A.colmajor[1][3] = 5;A.colmajor[2][3] = 10;A.colmajor[3][3] = 1;
	dtype b[4];
	b[0] = 4;
	b[1] = 9;
	b[2] = 8;
	b[3] = 10;
	dtype x_expected[4];
	x_expected[0] = 0.726028384122182;
	x_expected[1] = 0.126936099704632;
	x_expected[2] = 0.0;
	x_expected[3] = 0.231431453065341;

	NNLS_Single_Input input = NNLS_Single_Input(&A,b,10);
	int iterations = nnls_single_cpu(input);

	// printvectors(input.x,x_expected,4);
	checkvector(input.x,x_expected,
		4,"nnls_single_cpu",1e-5);
	printf("(%d iterations)\n",iterations);

}

void test_nnls_multiple_cpu_singlerhs()
{
	DenseMatrix A = DenseMatrix(4,4);
	A.colmajor[0][0] = 8;A.colmajor[1][0] = 6;A.colmajor[2][0] = 10;A.colmajor[3][0] = 10;
	A.colmajor[0][1] = 9;A.colmajor[1][1] = 1;A.colmajor[2][1] = 10;A.colmajor[3][1] = 5;
	A.colmajor[0][2] = 1;A.colmajor[1][2] = 3;A.colmajor[2][2] = 2;A.colmajor[3][2] = 8;
	A.colmajor[0][3] = 9;A.colmajor[1][3] = 5;A.colmajor[2][3] = 10;A.colmajor[3][3] = 1;
	dtype b[4];
	b[0] = 4;
	b[1] = 9;
	b[2] = 8;
	b[3] = 10;
	dtype x_expected[4];
	x_expected[0] = 0.726028384122182;
	x_expected[1] = 0.126936099704632;
	x_expected[2] = 0.0;
	x_expected[3] = 0.231431453065341;
	dtype** CTb = new dtype*[1];
	CTb[0] = new dtype[4];
	matvecmult_transpose_cpu(A,b,CTb[0]);

	LowerTriangularMatrix CTC = LowerTriangularMatrix(4);
	matmult_ata_lowertriangular_cpu(CTC,A);
	NNLS_Multiple_Input input = NNLS_Multiple_Input(&CTC,CTb,1,7);

	int iterations = nnls_multiple_cpu(input);

	// printvectors(input.X[0],x_expected,4);
	checkvector(input.X[0],x_expected,
		4,"test_nnls_multiple_cpu_singlerhs",1e-5);
	printf("(%d iterations)\n",iterations);
	delete[] CTb[0];
	delete[] CTb;
}
void test_nnls_multiple_cpu()
{
	DenseMatrix A = DenseMatrix(4,3);
	A.colmajor[0][0] = 8;A.colmajor[1][0] = 6;A.colmajor[2][0] = 10;
	A.colmajor[0][1] = 9;A.colmajor[1][1] = 1;A.colmajor[2][1] = 10;
	A.colmajor[0][2] = 1;A.colmajor[1][2] = 3;A.colmajor[2][2] = 2;
	A.colmajor[0][3] = 9;A.colmajor[1][3] = 5;A.colmajor[2][3] = 10;
	DenseMatrix B = DenseMatrix(4,3);
	B.colmajor[0][0] = 10;B.colmajor[1][0] = 4;B.colmajor[2][0] = 7;
	B.colmajor[0][1] = 5;B.colmajor[1][1] = 9;B.colmajor[2][1] = 0;
	B.colmajor[0][2] = 8;B.colmajor[1][2] = 8;B.colmajor[2][2] = 8;
	B.colmajor[0][3] = 1;B.colmajor[1][3] = 10;B.colmajor[2][3] = 9;
	dtype**CTB;CTB = new dtype*[3];
	dtype **X_expected;X_expected = new dtype*[3];
	for(int col=0;col<3;++col)
	{
		CTB[col] = new dtype[4];
		X_expected[col] = new dtype[3];
		matvecmult_transpose_cpu(A,B.colmajor[col],CTB[col]);
	}
	X_expected[0][0] = 0;
	X_expected[0][1] = 1.121233356692361;
	X_expected[0][2] = 0.114225648213034;
	X_expected[1][0] = 0;
	X_expected[1][1] = 0.268395234758234;
	X_expected[1][2] = 0.697967764540995;
	X_expected[2][0] = 0;
	X_expected[2][1] = 1.563380281690141;
	X_expected[2][2] = 0;
	LowerTriangularMatrix CTC = LowerTriangularMatrix(3);
	matmult_ata_lowertriangular_cpu(CTC,A);
	NNLS_Multiple_Input input = NNLS_Multiple_Input(&CTC,CTB,3,10);

	int iterations = nnls_multiple_cpu(input);

	checkvector(input.X[0],X_expected[0],
		3,"test_nnls_multiple_cpu, 1 of 3",1e-5);
	checkvector(input.X[1],X_expected[1],
		3,"test_nnls_multiple_cpu, 2 of 3",1e-5);
	checkvector(input.X[2],X_expected[2],
		3,"test_nnls_multiple_cpu, 3 of 3",1e-5);
	printf("(%d iterations max)\n",iterations);
	for(int col=0;col<3;++col)
	{
		delete[] CTB[col];
		delete[] X_expected[col];
	}
	delete[] CTB;
	delete[] X_expected;
}



void test_nmf_cpu()
{
	DenseMatrix A = DenseMatrix(5,4);
	DenseMatrix W = DenseMatrix(5,2);
	DenseMatrix W_expected = DenseMatrix(5,2);
	DenseMatrix H = DenseMatrix(2,4);
	DenseMatrix H_expected = DenseMatrix(2,4);
	A.colmajor[0][0] = 8;A.colmajor[1][0] = 1;A.colmajor[2][0] = 2;A.colmajor[3][0] = 1;
	A.colmajor[0][1] = 9;A.colmajor[1][1] = 3;A.colmajor[2][1] = 10;A.colmajor[3][1] = 4;
	A.colmajor[0][2] = 1;A.colmajor[1][2] = 5;A.colmajor[2][2] = 10;A.colmajor[3][2] = 9;
	A.colmajor[0][3] = 9;A.colmajor[1][3] = 10;A.colmajor[2][3] = 5;A.colmajor[3][3] = 8;
	A.colmajor[0][4] = 6;A.colmajor[1][4] = 10;A.colmajor[2][4] = 8;A.colmajor[3][4] = 10;
	A.copyColumnToRow();
	// printcolmajor(A.colmajor,A.rows,A.cols);
	// printrowmajor(A.rowmajor,A.rows,A.cols);
       
    //Only approximate! Forgot to 'format long' in MATLAB when copying these numbers.
	W_expected.colmajor[0][0] = 0;W_expected.colmajor[1][0] = 8.3240;
	W_expected.colmajor[0][1] = 6.9608;W_expected.colmajor[1][1] = 9.9018;
	W_expected.colmajor[0][2] = 13.6865;W_expected.colmajor[1][2] = 1.1497;
	W_expected.colmajor[0][3] = 10.5892;W_expected.colmajor[1][3] = 9.2789;
	W_expected.colmajor[0][4] = 14.3970;W_expected.colmajor[1][4] = 6.2455;
      
	H_expected.colmajor[0][0] = 0;H_expected.colmajor[1][0] = 0.5120;
	H_expected.colmajor[0][1] = 0.9458;H_expected.colmajor[1][1] = 0.2017;
	H_expected.colmajor[2][0] = 0.5520;H_expected.colmajor[3][0] = 0.6581;
	H_expected.colmajor[2][1] = 0.2492;H_expected.colmajor[3][1] = 0.0527;


	NMF_Input input = NMF_Input(&W,&H,&A,100,10);
	
	nmf_cpu(input);

	// printcolmajor(W.colmajor,W.rows,W.cols);
	// printcolmajor(H.colmajor,H.rows,H.cols);

	// checkMatrix(W.colmajor,W_expected.colmajor,2,5,"nmf_cpu, W",1e-5);
	// checkMatrix(H.colmajor,H_expected.colmajor,4,2,"nmf_cpu, H",1e-5);

	W.copyColumnToRow();
	printf("Calculated solution approximate Frobenius norm: %f\n",Fnorm(W,H,A));
	W_expected.copyColumnToRow();
	printf("MATLAB nnmf() approximate Frobenius norm: %f\n",Fnorm(W_expected,H_expected,A));
	printf("test_nmf_cpu completed.\n");
}


void test_reader(std::string filename,const char delimiter)
{
	DenseMatrix* A = readMatrix(filename,delimiter);
	printf("Completed reader test.\n");
	if(A) delete A;
}


void test_matmatmult_colmajor_cpu()
{
	int Arows = 5;
	int Acols = 3;
	int Brows = Acols;
	int Bcols = 4;
	int Crows = Arows;
	int Ccols = Bcols;

	dtype**A_colmajor = new dtype*[Acols];
	for(int col = 0; col < Acols; ++col) A_colmajor[col] = new dtype[Arows];
	A_colmajor[0][0]=2;A_colmajor[1][0]=6;A_colmajor[2][0]=6;
	A_colmajor[0][1]=8;A_colmajor[1][1]=8;A_colmajor[2][1]=7;
	A_colmajor[0][2]=0;A_colmajor[1][2]=7;A_colmajor[2][2]=7;
	A_colmajor[0][3]=3;A_colmajor[1][3]=9;A_colmajor[2][3]=2;
	A_colmajor[0][4]=7;A_colmajor[1][4]=1;A_colmajor[2][4]=5;
	dtype**B_colmajor = new dtype*[Bcols];
	for(int col = 0; col < Bcols; ++col) B_colmajor[col] = new dtype[Brows];
	B_colmajor[0][0]=2;B_colmajor[1][0]=9;B_colmajor[2][0]=9;B_colmajor[3][0]=4;
	B_colmajor[0][1]=2;B_colmajor[1][1]=7;B_colmajor[2][1]=2;B_colmajor[3][1]=3;
	B_colmajor[0][2]=2;B_colmajor[1][2]=3;B_colmajor[2][2]=3;B_colmajor[3][2]=6;
	dtype**C_colmajor = new dtype*[Ccols];
	for(int col = 0; col < Ccols; ++col) C_colmajor[col] = new dtype[Crows];
	C_colmajor[0][0]=28;C_colmajor[1][0]=78;C_colmajor[2][0]=48;C_colmajor[3][0]=62;
	C_colmajor[0][1]=46;C_colmajor[1][1]=149;C_colmajor[2][1]=109;C_colmajor[3][1]=98;
	C_colmajor[0][2]=28;C_colmajor[1][2]=70;C_colmajor[2][2]=35;C_colmajor[3][2]=63;
	C_colmajor[0][3]=28;C_colmajor[1][3]=96;C_colmajor[2][3]=51;C_colmajor[3][3]=51;
	C_colmajor[0][4]=26;C_colmajor[1][4]=85;C_colmajor[2][4]=80;C_colmajor[3][4]=61;

	// randInit(A_colmajor,Arows,Acols);
	// randInit(B_colmajor,Brows,Bcols);

	DenseMatrix A = DenseMatrix(A_colmajor,Arows,Acols);
	DenseMatrix B = DenseMatrix(B_colmajor,Brows,Bcols);
	DenseMatrix C = DenseMatrix(Crows,Ccols);

	matmatmult_colmajor_cpu(A,B,C);
	
	printcolmajor(C_colmajor,C.rows,C.cols);
	std::cout << std::endl;
	printcolmajor(C.colmajor,C.rows,C.cols);
	checkMatrix(C.colmajor,C_colmajor,C.cols,C.rows,"test_matmatmult_colmajor_cpu",1e-5);

	for(int col = 0; col < Ccols; ++col) delete[] C_colmajor[col];
	delete[] C_colmajor;
}
