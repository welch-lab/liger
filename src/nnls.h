#pragma once

#include <math.h>
#include <ctime>
#include <vector>

#include "ls.h"
#include "mask.h"

void advanceInfeasibilityState_single(NNLS_Single_State& state)
{
	if( state.infeasible < state.lowest_infeasible ) 
	{
		state.lowest_infeasible = state.infeasible;
		state.full_exchange_buffer = FULL_EXCHANGE_BUFFER;
		state.full_exchange_mode = true;
	}
	else
	{
		--state.full_exchange_buffer;
		if ( state.full_exchange_buffer <= 0 )
		{
			state.full_exchange_mode = false;
		}
	}
}

int nnls_single_cpu(NNLS_Single_Input& input)
{//nonnegative least squares of Cx - b
	int cols = input.C->cols;
	int rows = input.C->rows;
	NNLS_Single_State state = NNLS_Single_State(rows,cols);

	matvecmult_transpose_cpu(*input.C,input.b,state.y_masked,-1);
	state.infeasible = markInfeasible(*state.infeasiblemask,state.x_masked,state.y_masked,*state.xmask);
	while( state.infeasible > 0 and state.iterations < input.max_iter )
	{
		advanceInfeasibilityState_single(state);

		switchSets(*state.infeasiblemask,*state.xmask,state.full_exchange_mode);
		// printf("infeasiblemask: %d,%d,%d,%d\n",state.infeasiblemask->value[0],state.infeasiblemask->value[1],
		// 	state.infeasiblemask->value[2],state.infeasiblemask->value[3]);
		applyColumnMask(*input.C,*state.C_xmask,*state.xmask);
		applyColumnMask(*input.C,*state.C_ymask,*state.xmask,true);
		
		normal_equations_cpu(*state.C_xmask,state.x_masked,input.b);
		// printf("x_masked, xmask: %f,%f,%f,%f; %d,%d,%d,%d\n",state.x_masked[0],state.x_masked[1],state.x_masked[2],state.x_masked[3],
		// 	state.xmask->value[0],state.xmask->value[1],state.xmask->value[2],state.xmask->value[3]);

		matvecmult_colmajor_cpu(*state.C_xmask,state.x_masked,state.y_masked_intermediate);
		vectorsub(input.b,state.y_masked_intermediate,rows);
		matvecmult_transpose_cpu(*state.C_ymask,state.y_masked_intermediate,state.y_masked);
		// printf("y_masked: %f,%f,%f,%f\n",state.y_masked[0],state.y_masked[1],state.y_masked[2],state.y_masked[3]);
		
		state.infeasible = markInfeasible(*state.infeasiblemask,state.x_masked,state.y_masked,*state.xmask);
		++state.iterations;
	}

	overwriteOriginalWithMask(input.x,state.x_masked,*state.xmask);

	return state.iterations;
}

void initialize_multiple_cpu(NNLS_Multiple_Input& input,NNLS_Multiple_State& state)
{
	for(int i=0;i<state.cols_rhs;++i)
	{
		vectorinit(state.x_masked[i],state.cols,0);
		vectornegate(input.CTB[i],state.y_masked[i],state.cols);
		state.full_exchange_buffer[i] = FULL_EXCHANGE_BUFFER;
		state.lowest_infeasible[i] = state.cols + 1;
		state.full_exchange_mode[i] = true;
		state.infeasible[i] = state.cols+1;
		for(int j=0;j<state.xmasks[i]->dim;++j) state.xmasks[i]->value[j] = false;
		for(int j=0;j<state.infeasiblemasks[i]->dim;++j) state.infeasiblemasks[i]->value[j] = false;
		state.CFTCF[i]->dim = state.CGTCF[i]->originalcols;
		state.CFTCF[i]->totalsize = (state.CGTCF[i]->originalcols*(state.CGTCF[i]->originalcols+1))/2;
		state.CGTCF[i]->rows = state.CGTCF[i]->originalrows;
		state.CGTCF[i]->cols = state.CGTCF[i]->originalcols;
		state.CGTCF[i]->totalsize = state.CGTCF[i]->originalrows*state.CGTCF[i]->originalcols;
		//printf("%f,%f,%f,%f\n",state.y_masked[i][0],state.y_masked[i][1],state.y_masked[i][2],state.y_masked[i][3]);
	}
	state.totalfeasible = 0;
	state.iterations = 0;
	destroyCholeskyMap(state.choleskyMap);
}
void markInfeasible_multiple_cpu(NNLS_Multiple_State& state)
{
	for(int i=0;i<state.cols_rhs;++i)
	{
		if(state.infeasible[i] > 0)
		{
			state.infeasible[i] = markInfeasible(*state.infeasiblemasks[i],state.x_masked[i],state.y_masked[i],*state.xmasks[i]);
			if(state.infeasible[i] == 0) ++state.totalfeasible;
		}
	}
}

void determineCholeskyFactors_cpu(NNLS_Multiple_Input& input,NNLS_Multiple_State& state)
{
	for(int i=0;i<state.cols_rhs;++i)
	{
		if(state.infeasible[i] > 0)
		{
			std::string key = maskToString(*state.xmasks[i]);
			CholeskyMap::iterator it = state.choleskyMap.find(key);
			if( it != state.choleskyMap.end() )
			{
				state.G[i] = it->second;
			}
			else
			{
				generateCFTCF(*state.CFTCF[i],*input.CTC,*state.xmasks[i]);
				LowerTriangularMatrix* G = new LowerTriangularMatrix(state.CFTCF[i]->dim);
				// printf("state.CFTCF[i]->dim: %d\n",state.CFTCF[i]->dim);
				// printf("G dim,totalsize: %d,%d\n",G->dim,G->totalsize);
				// std::cout << "G key: " << key << std::endl;
				cholesky_lowertriangular_cpu(*G,*state.CFTCF[i]);
				state.G[i] = G;
				std::pair<CholeskyMap::iterator,bool> ret = 
					state.choleskyMap.insert ( std::pair<std::string,LowerTriangularMatrix*>(key,G) );
				if ( not ret.second ) std::cout << "ERROR! Duplicate Cholesky factor was inserted." << std::endl;
			}
		}
	}
}
void determineCholeskyFactors_cpu_profile(NNLS_Multiple_Input& input,NNLS_Multiple_State& state)
{
	std::clock_t start;
	std::clock_t checkStart = std::clock();
	// #pragma omp parallel for
	for(int i=0;i<state.cols_rhs;++i)
	{
		if(state.infeasible[i] > 0)
		{
			start = std::clock();
			// maskToString3(*state.xmasks[i]);
			// std::string& key = state.xmasks[i]->tmp;
			// maskToString4(*state.xmasks[i]);std::string& key = state.xmasks[i]->tmp;
			maskToString5(*state.xmasks[i]);std::string key = &(*(state.xmasks[i]->vs.begin()));
			// std::string key; maskToString2(*state.xmasks[i],key);
			// std::vector<int> key; maskToVector(*state.xmasks[i],key);
			input.choleskyMaskToString_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
			
			start = std::clock();
			CholeskyMap::iterator it = state.choleskyMap.find(key);
			// CholeskyMap2::iterator it = state.choleskyMap2.find(key);
			input.choleskyKeyFind_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

			// if( it != state.choleskyMap2.end() ) state.G[i] = it->second;
			if( it != state.choleskyMap.end() ) state.G[i] = it->second;
			else
			{
				start = std::clock();
				generateCFTCF(*state.CFTCF[i],*input.CTC,*state.xmasks[i]);
				input.choleskyCFTCF_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

				start = std::clock();
				LowerTriangularMatrix* G = new LowerTriangularMatrix(state.CFTCF[i]->dim);
				input.choleskyAllocate_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

				// printf("state.CFTCF[i]->dim: %d\n",state.CFTCF[i]->dim);
				// printf("G dim,totalsize: %d,%d\n",G->dim,G->totalsize);
				// std::cout << "G key: " << key << std::endl;
				
				start = std::clock();
				cholesky_lowertriangular_cpu(*G,*state.CFTCF[i]);
				input.choleskyFactorization_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

				state.G[i] = G;
				start = std::clock();
				std::pair<CholeskyMap::iterator,bool> ret = 
					state.choleskyMap.insert ( std::pair<std::string,LowerTriangularMatrix*>(key,G) );
				// std::pair<CholeskyMap2::iterator,bool> ret = 
				// 	state.choleskyMap2.insert ( std::pair<std::vector<int>,LowerTriangularMatrix*>(key,G) );
				if ( not ret.second ) std::cout << "ERROR! Duplicate Cholesky factor was inserted." << std::endl;
				input.choleskyMapInsert_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
			}
		}
	}
	input.choleskyCheck_time += ( std::clock() - checkStart ) / (double) CLOCKS_PER_SEC;
}

void advanceInfeasibilityState_multiple(NNLS_Multiple_State& state)
{
	for(int i=0;i<state.cols_rhs;++i)
	{
		if( state.infeasible[i] < state.lowest_infeasible[i] ) 
		{
			state.lowest_infeasible[i] = state.infeasible[i];
			state.full_exchange_buffer[i] = FULL_EXCHANGE_BUFFER;
			state.full_exchange_mode[i] = true;
		}
		else
		{
			--state.full_exchange_buffer[i];
			if ( state.full_exchange_buffer[i] <= 0 )
			{
				state.full_exchange_mode[i] = false;
			}
		}
	}
}

int nnls_multiple_cpu(NNLS_Multiple_Input& input)
{//nonnegative least squares of CX - B
	// int cols = input.CTC->dim;
	// int cols_rhs = input.cols_rhs;
	// NNLS_Multiple_State state = NNLS_Multiple_State(cols,cols_rhs);
	NNLS_Multiple_State& state = *input.state;

	initialize_multiple_cpu(input,state);

	markInfeasible_multiple_cpu(state);
	determineCholeskyFactors_cpu(input,state);

	while( state.totalfeasible < state.cols_rhs and state.iterations < input.max_iter )
	{
		advanceInfeasibilityState_multiple(state);

		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) switchSets(*state.infeasiblemasks[i],*state.xmasks[i],state.full_exchange_mode[i]);
		// printf("infeasiblemask: %d,%d,%d,%d\n",state.infeasiblemasks[0]->value[0],state.infeasiblemasks[0]->value[1],
		// 	state.infeasiblemasks[0]->value[2],state.infeasiblemasks[0]->value[3]);
		determineCholeskyFactors_cpu(input,state);
		
		// printv(state.G[0]->rowmajor,state.G[0]->totalsize);
		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) applyVectorMask(input.CTB[i],state.x_masked[i],*state.xmasks[i]);
		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) normal_equations_precomputedCholesky_cpu(*state.G[i],state.x_masked[i]);
		// printf("x_masked, xmask: %f,%f,%f,%f; %d,%d,%d,%d\n",state.x_masked[0][0],state.x_masked[0][1],state.x_masked[0][2],state.x_masked[0][3],
		// 	state.xmasks[0]->value[0],state.xmasks[0]->value[1],state.xmasks[0]->value[2],state.xmasks[0]->value[3]);

		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) generateCGTCF(*state.CGTCF[i],*input.CTC,*state.xmasks[i]);
		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) matvecmult_colmajor_cpu(*state.CGTCF[i],state.x_masked[i],state.y_masked[i]);
		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) generateCGTb(input.CTB[i],state.CGTb[i],*state.xmasks[i]);
		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) vectorsub(state.CGTb[i],state.y_masked[i],state.CGTCF[i]->rows);
		// printf("y_masked: %f,%f,%f,%f\n",state.y_masked[0][0],state.y_masked[0][1],state.y_masked[0][2],state.y_masked[0][3]);
		
		markInfeasible_multiple_cpu(state);
		++state.iterations;
	}

	for(int i=0;i<state.cols_rhs;++i) overwriteOriginalWithMask(input.X[i],state.x_masked[i],*state.xmasks[i]);

	int totaliterations = state.iterations;

	return totaliterations;
}

int nnls_multiple_cpu_profile(NNLS_Multiple_Input& input)
{//nonnegative least squares of CX - B
	std::clock_t start;

	// start = std::clock();
	// int cols = input.CTC->dim;
	// int cols_rhs = input.cols_rhs;
	// NNLS_Multiple_State state = NNLS_Multiple_State(cols,cols_rhs);
	// input.allocate_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	NNLS_Multiple_State& state = *input.state;

	start = std::clock();
	initialize_multiple_cpu(input,state);
	input.init_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

	start = std::clock();
	markInfeasible_multiple_cpu(state);
	input.mark_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	
	start = std::clock();
	determineCholeskyFactors_cpu_profile(input,state);
	input.determine_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

	while( state.totalfeasible < state.cols_rhs and state.iterations < input.max_iter )
	{
		start = std::clock();
		advanceInfeasibilityState_multiple(state);
		input.advance_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

		start = std::clock();
		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) switchSets(*state.infeasiblemasks[i],*state.xmasks[i],state.full_exchange_mode[i]);
		input.switch_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		// printf("infeasiblemask: %d,%d,%d,%d\n",state.infeasiblemasks[0]->value[0],state.infeasiblemasks[0]->value[1],
		// 	state.infeasiblemasks[0]->value[2],state.infeasiblemasks[0]->value[3]);
		start = std::clock();
		determineCholeskyFactors_cpu_profile(input,state);
		input.determine_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		
		// printv(state.G[0]->rowmajor,state.G[0]->totalsize);
		start = std::clock();
		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) applyVectorMask(input.CTB[i],state.x_masked[i],*state.xmasks[i]);
		input.apply_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

		start = std::clock();
		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) normal_equations_precomputedCholesky_cpu(*state.G[i],state.x_masked[i]);
		input.normal_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		// printf("x_masked, xmask: %f,%f,%f,%f; %d,%d,%d,%d\n",state.x_masked[0][0],state.x_masked[0][1],state.x_masked[0][2],state.x_masked[0][3],
		// 	state.xmasks[0]->value[0],state.xmasks[0]->value[1],state.xmasks[0]->value[2],state.xmasks[0]->value[3]);

		start = std::clock();
		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) generateCGTCF(*state.CGTCF[i],*input.CTC,*state.xmasks[i]);
		input.generateCGTCF_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		
		start = std::clock();
		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) matvecmult_colmajor_cpu(*state.CGTCF[i],state.x_masked[i],state.y_masked[i]);
		input.matvec_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		
		start = std::clock();
		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) generateCGTb(input.CTB[i],state.CGTb[i],*state.xmasks[i]);
		input.generateCGTb_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

		start = std::clock();
		for(int i=0;i<state.cols_rhs;++i) if(state.infeasible[i] > 0) vectorsub(state.CGTb[i],state.y_masked[i],state.CGTCF[i]->rows);
		input.vector_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		// printf("y_masked: %f,%f,%f,%f\n",state.y_masked[0][0],state.y_masked[0][1],state.y_masked[0][2],state.y_masked[0][3]);

		start = std::clock();		
		markInfeasible_multiple_cpu(state);
		input.mark_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

		++state.iterations;
	}

	start = std::clock();
	for(int i=0;i<state.cols_rhs;++i) overwriteOriginalWithMask(input.X[i],state.x_masked[i],*state.xmasks[i]);
	input.overwrite_time += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

	int totaliterations = state.iterations;

	return totaliterations;
}