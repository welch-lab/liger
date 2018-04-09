#pragma once

#include "entities.h"
#include <stdio.h>
#include <stdlib.h>

typedef struct gpu_input
{
	dtype *d_W_colmajor,*d_H_colmajor,*d_A_colmajor;
	dtype *d_WTA_colmajor,*d_HAT_colmajor;
} gpu_input;

void nmf_gpu_profile(NMF_Input& input,gpu_input& gpu_in);

// void initHandle(gpu_input& in);
// void destroyHandle(gpu_input& in);
void colmatrix2gpu(dtype*device_colmatrix,dtype**host_colmatrix,int cols,int ld);
void colmatrix2host(dtype**host_colmatrix,dtype*device_colmatrix,int cols,int ld);
void gpufree(dtype*d_vec);
void gpumalloc(dtype*&d_vec,int size);
void gpusync();
void gpucheck(std::string msg);
// #define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
// inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
// {
//    if (code != cudaSuccess) 
//    {
//       fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
//       if (abort) exit(code);
//    }
// }

void level3_float(const float*d_A_colmajor,const float*d_B_colmajor,float*d_C_colmajor,
	int A_rows,int A_cols,int B_rows,int B_cols);
void level3_HAT_float(const float*d_H_colmajor,const float*d_A_colmajor,float*d_HAT_colmajor,
	int H_rows,int H_cols,int A_rows,int A_cols);
void level3_HAT_double(const double*d_H_colmajor,const double*d_A_colmajor,double*d_HAT_colmajor,
	int H_rows,int H_cols,int A_rows,int A_cols);
void level3_WTA_float(const float*d_W_colmajor,const float*d_A_colmajor,float*d_WTA_colmajor,
	int W_rows,int W_cols,int A_rows,int A_cols);
void level3_WTA_double(const double*d_W_colmajor,const double*d_A_colmajor,double*d_WTA_colmajor,
	int W_rows,int W_cols,int A_rows,int A_cols);