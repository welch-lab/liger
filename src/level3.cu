#include "parameters.h"
#include "entities.h"
#include "cublas_v2.h"
#include <iostream>
//#include <typeinfo>

#include "common.h"


// void initHandle(gpu_input& in)
// {
// 	cublasStatus_t status =	cublasCreate((cublasHandle_t*)in.handle);
// 	if( status != CUBLAS_STATUS_SUCCESS ) std::cout << "The cuBLAS handle was not properly created!\n" << std::endl;
// }
// void destroyHandle(gpu_input& in)
// {
// 	cublasStatus_t status =	cublasDestroy(*((cublasHandle_t*)in.handle));
// 	if( status != CUBLAS_STATUS_SUCCESS ) std::cout << "The cuBLAS handle was not properly destroyed!\n" << std::endl;
// }

void gpumalloc(dtype*& d_vec,int size)
{
	cudaError_t error;
	error = cudaMalloc((void **)& d_vec, size*sizeof(dtype));
	if( error != cudaSuccess ) std::cout << "GPU allocation failed!\n" << std::endl;
	// cudaDeviceSynchronize();
}
void vector2gpu(dtype*device_vector,dtype*host_vector,int size)
{
	cudaError_t error = 
		cudaMemcpy(device_vector,host_vector,size*sizeof(dtype),cudaMemcpyHostToDevice);
	if( error != cudaSuccess ) std::cout << "Host-to-device transfer failed! (vector2gpu) " << error << std::endl;
}
void colmatrix2gpu(dtype*device_colmatrix,dtype**host_colmatrix,int cols,int ld)
{
	for(int col=0;col<cols;++col)
	{
		int shift = col*ld;
		vector2gpu(device_colmatrix+shift,host_colmatrix[col],ld);
	}
}
void vector2host(dtype*host_vector,dtype*device_vector,int size)
{
	cudaError_t error = 
		cudaMemcpy(host_vector,device_vector,size*sizeof(dtype),cudaMemcpyDeviceToHost);
	if( error != cudaSuccess ) std::cout << "Device-to-host transfer failed! (vector2host) " << error << " " << size << std::endl;
}
void colmatrix2host(dtype**host_colmatrix,dtype*device_colmatrix,int cols,int ld)
{
	for(int col=0;col<cols;++col)
	{
		int shift = col*ld;
		vector2host(host_colmatrix[col],device_colmatrix+shift,ld);
		// std::cout << "shift, cols, ld: " << shift << ", " << cols << ", " << ld << std::endl;
	}
}

void gpufree(dtype*d_vec)
{
	cudaError_t error;
	error = cudaFree(d_vec);
	if( error != cudaSuccess ) std::cout << "cudaFree(W) failed!\n" << std::endl;
}
void gpusync()
{
	cudaDeviceSynchronize();
}
void gpucheck(std::string msg)
{
	cudaError_t error = cudaPeekAtLastError();
	std::cout << "Cuda Peek Error (" << msg << "): " << error << std::endl;
	std::cout << "cudaSuccess: " << cudaSuccess << std::endl;
}	







void level3_float(const float*d_A_colmajor,const float*d_B_colmajor,float*d_C_colmajor,
	int A_rows,int A_cols,int B_rows,int B_cols)
{
	cublasHandle_t handle;
	cublasStatus_t status = cublasCreate ( & handle );

	// if( H_cols != A_rows ) printf("ERROR! H_cols and A_rows are not equal!\n");
	float alpha = 1.0;
	float beta = 0.0;
	cublasSgemm(handle,
		CUBLAS_OP_N,CUBLAS_OP_N,
		A_rows,B_cols,A_cols,
		&alpha,
		d_A_colmajor,A_rows,
		d_B_colmajor,B_rows,
		&beta,
		d_C_colmajor,A_rows
		);

	cublasDestroy(handle);
}
void level3_HAT_float(const float*d_H_colmajor,const float*d_A_colmajor,float*d_HAT_colmajor,
	int H_rows,int H_cols,int A_rows,int A_cols)
{
	cublasHandle_t handle;
	cublasStatus_t status = cublasCreate ( & handle );

	// if( H_cols != A_rows ) printf("ERROR! H_cols and A_rows are not equal!\n");
	float alpha = 1.0;
	float beta = 0.0;
	cublasSgemm(handle,
		CUBLAS_OP_N,CUBLAS_OP_T,
		H_rows,A_rows,H_cols,
		&alpha,
		d_H_colmajor,H_rows,
		d_A_colmajor,A_rows,
		&beta,
		d_HAT_colmajor,H_rows
		);

	cublasDestroy(handle);
}
void level3_HAT_double(const double*d_H_colmajor,const double*d_A_colmajor,double*d_HAT_colmajor,
	int H_rows,int H_cols,int A_rows,int A_cols)
{
	cublasHandle_t handle;
	cublasStatus_t status = cublasCreate ( & handle );

	// if( H_cols != A_rows ) printf("ERROR! H_cols and A_rows are not equal!\n");
	double alpha = 1.0;
	double beta = 0.0;
	cublasDgemm(handle,
		CUBLAS_OP_N,CUBLAS_OP_T,
		H_rows,A_rows,H_cols,
		&alpha,
		d_H_colmajor,H_rows,
		d_A_colmajor,A_rows,
		&beta,
		d_HAT_colmajor,H_rows
		);

	cublasDestroy(handle);
}
void level3_WTA_float(const float*d_W_colmajor,const float*d_A_colmajor,float*d_WTA_colmajor,
	int W_rows,int W_cols,int A_rows,int A_cols)
{
	cublasHandle_t handle;
	cublasStatus_t status =	cublasCreate ( & handle );

	// if( W_cols != A_rows ) printf("ERROR! W_cols and A_rows are not equal!\n");
	float alpha = 1.0;
	float beta = 0.0;
	cublasSgemm(handle,
		CUBLAS_OP_T,CUBLAS_OP_N,
		W_cols,A_cols,W_rows,
		&alpha,
		d_W_colmajor,W_rows,
		d_A_colmajor,A_rows,
		&beta,
		d_WTA_colmajor,W_cols
		);

	cublasDestroy(handle);
}
void level3_WTA_double(const double*d_W_colmajor,const double*d_A_colmajor,double*d_WTA_colmajor,
	int W_rows,int W_cols,int A_rows,int A_cols)
{
	cublasHandle_t handle;
	cublasStatus_t status =	cublasCreate ( & handle );

	// if( W_cols != A_rows ) printf("ERROR! W_cols and A_rows are not equal!\n");
	double alpha = 1.0;
	double beta = 0.0;
	cublasDgemm(handle,
		CUBLAS_OP_T,CUBLAS_OP_N,
		W_cols,A_cols,W_rows,
		&alpha,
		d_W_colmajor,W_rows,
		d_A_colmajor,A_rows,
		&beta,
		d_WTA_colmajor,W_cols
		);

	cublasDestroy(handle);
}
