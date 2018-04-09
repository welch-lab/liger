#pragma once

#include "parameters.h"
#include "entities.h"


dtype averageMatrix(DenseMatrix* A)
{
	dtype sum = 0;
	for(int row = 0; row < A->rows; ++row)
	{
		for(int col = 0; col < A->cols; ++col)
		{
			sum += A->colmajor[col][row];
		}
	}
	return sum / ( A->rows*A->cols );
}
void initializeMatrix(DenseMatrix* M, dtype factor=1)
{
	// dtype num = 1;
	for(int col = 0;col < M->cols;++col)
	{
		for(int row=0;row<M->rows;++row)
		{
			dtype value = (rand()%2)*factor;
			M->colmajor[col][row] = value;
			M->rowmajor[row][col] = value;
		}
	}
}

void printr(dtype** M,int rows,int cols)
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
void printc(dtype** M,int rows,int cols)
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
void printl(LowerTriangularMatrix* L)
{
	for(int row=0;row<L->dim;++row)
	{
		for(int col=0;col<=row;++col)
		{
			printf("%f,",L->rowmajor[(row*(row+1))/2+col]);
		}
		printf("\n");
	}
}
void nmf_cpu(NMF_Input& input)
{
	NMF_State state = NMF_State(input.m,input.k,input.n);
	NNLS_Multiple_Input nnls_input_1 = NNLS_Multiple_Input(state.HHT,input.W->rowmajor,state.HAT,input.m,input.max_iter_nnls);
	NNLS_Multiple_Input nnls_input_2 = NNLS_Multiple_Input(state.WTW,input.H->colmajor,state.WTA,input.n,input.max_iter_nnls);

	// initializeMatrix(input.H,averageMatrix(input.A));
	initializeMatrix(input.H);
	int iterations = 0;
	int iterations_nnls = 0;
	while(iterations < input.max_iter_nmf)//TODO: implement solution-sensitive stopping criterion
	{
		input.H->copyColumnToRow();
		matmult_ata_lowertriangular_pointers_cpu(*state.HHT,input.H->rowmajor,input.H->cols);
		for(int i=0;i<input.m;++i) matvecmult_colmajor_cpu(*input.H,input.A->rowmajor[i],state.HAT[i]);

		// printf("H: \n");printr(input.H->rowmajor,input.H->rows,input.H->cols);
		// printc(state.HAT,input.k,input.m);
		// printl(state.HHT);
		// printf("HHT: \n");printl(nnls_input_1.CTC);
		// printf("HAT: \n");printc(nnls_input_1.CTB,input.k,input.m);
		// int iterations1 = 
		iterations_nnls += nnls_multiple_cpu(nnls_input_1);
		// printf("nnls iter 1: %d\n",iterations1);
		// printc(nnls_input_1.X,input.k,input.m);

		input.W->copyRowToColumn();
		matmult_ata_lowertriangular_cpu(*state.WTW,*input.W);
		for(int i=0;i<input.n;++i) matvecmult_transpose_cpu(*input.W,input.A->colmajor[i],state.WTA[i]);
		
		// printf("WTW: \n");printl(nnls_input_2.CTC);
		// printf("WTA: \n");printc(nnls_input_2.CTB,input.k,input.n);
		// int iterations2 = 
		iterations_nnls += nnls_multiple_cpu(nnls_input_2);
		// printf("nnls iter 2: %d\n",iterations2);
		// printc(nnls_input_2.X,input.k,input.n);

		++iterations;
	}

	printf("Average NNLS Iterations: %f\n",(dtype)iterations_nnls/(dtype)iterations/2.0);
	printf("NMF Iterations: %d\n",iterations);
}

#include <ctime>
void print_nnls_time(NNLS_Multiple_Input& in)
{
	double total_time = in.allocate_time+in.init_time+in.advance_time+in.switch_time
		+in.determine_time+in.apply_time+in.normal_time+in.generateCGTCF_time+in.matvec_time
		+in.generateCGTb_time+in.vector_time+in.mark_time+in.overwrite_time;
	printf("NNLS Profiling results (in seconds): (total: %f)\n",total_time);
	printf("\tAllocation time: %f (%f %%)\n",in.allocate_time,in.allocate_time*100.0/total_time);
	printf("\tInitialization time: %f (%f %%)\n",in.init_time,in.init_time*100.0/total_time);
	printf("\tAdvance infeasible state time: %f (%f %%)\n",in.advance_time,in.advance_time*100.0/total_time);
	printf("\tSwitch sets time: %f (%f %%)\n",in.switch_time,in.switch_time*100.0/total_time);
	printf("\tDetermine Cholesky factors time: %f (%f %%)\n",in.determine_time,in.determine_time*100.0/total_time);
		printf("\t\tCholesky mask to string time: %f (%f %% of total Cholesky time)\n",in.choleskyMaskToString_time,in.choleskyMaskToString_time*100.0/in.determine_time);
		printf("\t\tCholesky key find time: %f (%f %% of total Cholesky time)\n",in.choleskyKeyFind_time,in.choleskyKeyFind_time*100.0/in.determine_time);
		printf("\t\tCholesky CFTCF time: %f (%f %% of total Cholesky time)\n",in.choleskyCFTCF_time,in.choleskyCFTCF_time*100.0/in.determine_time);
		printf("\t\tCholesky allocate time: %f (%f %% of total Cholesky time)\n",in.choleskyAllocate_time,in.choleskyAllocate_time*100.0/in.determine_time);
		printf("\t\tCholesky factorization time: %f (%f %% of total Cholesky time)\n",in.choleskyFactorization_time,in.choleskyFactorization_time*100.0/in.determine_time);
		printf("\t\tCholesky map insert time: %f (%f %% of total Cholesky time)\n",in.choleskyMapInsert_time,in.choleskyMapInsert_time*100.0/in.determine_time);
		printf("\t\tCholesky check time: %f (%f %% of total Cholesky time)\n",in.choleskyCheck_time,in.choleskyCheck_time*100.0/in.determine_time);
	printf("\tApply masks time: %f (%f %%)\n",in.apply_time,in.apply_time*100.0/total_time);
	printf("\tNormal equations least squares time: %f (%f %%)\n",in.normal_time,in.normal_time*100.0/total_time);
	printf("\tGenerate CGTCF time: %f (%f %%)\n",in.generateCGTCF_time,in.generateCGTCF_time*100.0/total_time);
	printf("\tCGCTF*x matvec time: %f (%f %%)\n",in.matvec_time,in.matvec_time*100.0/total_time);
	printf("\tGenerate CGTb time: %f (%f %%)\n",in.generateCGTb_time,in.generateCGTb_time*100.0/total_time);
	printf("\tCGTCF*x-CGTb vectorsub time: %f (%f %%)\n",in.vector_time,in.vector_time*100.0/total_time);
	printf("\tMark infeasible time: %f (%f %%)\n",in.mark_time,in.mark_time*100.0/total_time);
	printf("\tSolution overwrite time: %f (%f %%)\n",in.overwrite_time,in.overwrite_time*100.0/total_time);
}
void nmf_cpu_profile(NMF_Input& input)
{
	std::clock_t start;
	double allocation_time=0,H_init_time=0,H_copy_time = 0,W_copy_time = 0,HHT_time=0,WTW_time=0,HAT_time=0,WTA_time=0,WT_solve_time=0,H_solve_time=0;

	start = std::clock();
	NMF_State state = NMF_State(input.m,input.k,input.n);
	NNLS_Multiple_Input nnls_input_1 = NNLS_Multiple_Input(state.HHT,input.W->rowmajor,state.HAT,input.m,input.max_iter_nnls);
	NNLS_Multiple_Input nnls_input_2 = NNLS_Multiple_Input(state.WTW,input.H->colmajor,state.WTA,input.n,input.max_iter_nnls);
	// NNLS_Multiple_State nnls_state_1 = NNLS_Multiple_State(nnls_input_1.CTC->dim,nnls_input_1.cols_rhs);
	// NNLS_Multiple_State nnls_state_2 = NNLS_Multiple_State(nnls_input_2.CTC->dim,nnls_input_2.cols_rhs);
	// nnls_input_1.state = &nnls_state_1;
	// nnls_input_2.state = &nnls_state_2;
	allocation_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

	// initializeMatrix(input.H,averageMatrix(input.A));
	start = std::clock();
	initializeMatrix(input.H);
	H_init_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

	int iterations = 0;
	int iterations_nnls = 0;
	while(iterations < input.max_iter_nmf)//TODO: implement solution-sensitive stopping criterion
	{
		start = std::clock();
		input.H->copyColumnToRow();
		H_copy_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

		start = std::clock();
		matmult_ata_lowertriangular_pointers_cpu(*state.HHT,input.H->rowmajor,input.H->cols);
		HHT_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

		start = std::clock();
		for(int i=0;i<input.m;++i) matvecmult_colmajor_cpu(*input.H,input.A->rowmajor[i],state.HAT[i]);
		HAT_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

		// printf("H: \n");printr(input.H->rowmajor,input.H->rows,input.H->cols);
		// printc(state.HAT,input.k,input.m);
		// printl(state.HHT);
		// printf("HHT: \n");printl(nnls_input_1.CTC);
		// printf("HAT: \n");printc(nnls_input_1.CTB,input.k,input.m);
		// int iterations1 = 
		start = std::clock();
		iterations_nnls += nnls_multiple_cpu_profile(nnls_input_1);
		WT_solve_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		// printf("nnls iter 1: %d\n",iterations1);
		// printc(nnls_input_1.X,input.k,input.m);

		start = std::clock();
		input.W->copyRowToColumn();
		W_copy_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

		start = std::clock();
		matmult_ata_lowertriangular_cpu(*state.WTW,*input.W);
		WTW_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

		start = std::clock();
		for(int i=0;i<input.n;++i) matvecmult_transpose_cpu(*input.W,input.A->colmajor[i],state.WTA[i]);
		WTA_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		
		// printf("WTW: \n");printl(nnls_input_2.CTC);
		// printf("WTA: \n");printc(nnls_input_2.CTB,input.k,input.n);
		// int iterations2 = 
		start = std::clock();
		iterations_nnls += nnls_multiple_cpu_profile(nnls_input_2);
		H_solve_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		// printf("nnls iter 2: %d\n",iterations2);
		// printc(nnls_input_2.X,input.k,input.n);

		++iterations;
	}
	
	printf("\tAverage NNLS Iterations: %f\n",(dtype)iterations_nnls/(dtype)iterations/2.0);
	printf("\tNMF Iterations: %d\n",iterations);

	double total_time = allocation_time+H_init_time+
		H_copy_time+W_copy_time+
		H_solve_time+WT_solve_time+
		HAT_time+WTA_time+
		HHT_time+WTW_time;
	printf("NMF Profiling results (in seconds): (total: %f)\n",total_time);
	printf("\tState and input allocation time: %f (%f %%)\n",allocation_time,allocation_time*100.0/total_time);
	printf("\tH initialization time: %f (%f %%)\n",H_init_time,H_init_time*100.0/total_time);
	printf("\tH, W copy times: %f (%f %%), %f (%f %%)\n",H_copy_time,H_copy_time*100.0/total_time,W_copy_time,W_copy_time*100.0/total_time);
	printf("\tH, WT solve times (nnls): %f (%f %%), %f (%f %%)\n",H_solve_time,H_solve_time*100.0/total_time,WT_solve_time,WT_solve_time*100.0/total_time);
	printf("\tHAT, WTA times (matmult): %f (%f %%), %f (%f %%)\n",HAT_time,HAT_time*100.0/total_time,WTA_time,WTA_time*100.0/total_time);
	printf("\tHHT, WTW times (matmult): %f (%f %%), %f (%f %%)\n",HHT_time,HHT_time*100.0/total_time,WTW_time,WTW_time*100.0/total_time);
	print_nnls_time(nnls_input_1);
	print_nnls_time(nnls_input_2);

}