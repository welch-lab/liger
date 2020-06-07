#pragma once

#include <string>
#include <vector>

#include <map>
#if __cplusplus > 199711L
	#include <unordered_map>
#endif

#include "parameters.h" 

using std::endl;

typedef struct DenseMatrix
{
	dtype **rowmajor,**colmajor;
	int rows,cols,totalsize;
	bool dependent;
	int originalrows,originalcols;
	DenseMatrix();
	//DenseMatrix() : rows(0),cols(0),totalsize(0),dependent(false),rowmajor(NULL),colmajor(NULL) {}
	
	DenseMatrix(int rows_,int cols_,bool dependent_=false) : 
		rows(rows_),cols(cols_),totalsize(rows*cols),dependent(dependent_),
		originalrows(rows_),originalcols(cols_)
	{
		rowmajor = new dtype*[rows];
		colmajor = new dtype*[cols];
		if (not dependent)
		{
			for(int row=0;row<rows;++row) { rowmajor[row] = new dtype[cols](); }
			for(int col=0;col<cols;++col) { colmajor[col] = new dtype[rows](); }
		}
	}
	DenseMatrix(dtype**colmajor_,int rows_,int cols_,bool dependent_=false) : 
		colmajor(colmajor_),rows(rows_),cols(cols_),totalsize(rows*cols),dependent(dependent_),
		originalrows(rows_),originalcols(cols_)
	{
		rowmajor = new dtype*[rows];
		if (not dependent)
		{
			for(int row=0;row<rows;++row) { rowmajor[row] = new dtype[cols](); }
		}
	}
	void copyColumnToRow()
	{
		for(int col=0;col<cols;++col)
		{
			for(int row=0;row<rows;++row)
			{
				rowmajor[row][col] = colmajor[col][row];
			}
		}
	}
	void copyRowToColumn()
	{
		for(int row=0;row<rows;++row)
		{
			for(int col=0;col<cols;++col)
			{
				colmajor[col][row] = rowmajor[row][col];
			}
		}
	}
	void resetDim(int rows_,int cols_)
	{
		rows = rows_;
		cols = cols_;
		totalsize = rows*cols;
	}
	~DenseMatrix()
	{
		if(not dependent)
		{
			for(int row=0;row<originalrows;++row) { delete[] rowmajor[row]; }
			for(int col=0;col<originalcols;++col) { delete[] colmajor[col]; }
		}
		delete[] rowmajor;
		delete[] colmajor;
	}
} DenseMatrix;

DenseMatrix operator +(const DenseMatrix& X, const DenseMatrix& Y) {
	if (X.rows != Y.rows || X.cols != Y.cols)
	{
		//mismatched dimensions in matrix addition. Returning X
		return X;
	}
	DenseMatrix Z = DenseMatrix(X.rows,X.cols,true);
	for(int i=0; i < Z.rows; ++i)
	{
		for(int j=0; j < Z.cols; ++j)
		{
			Z.rowmajor[i][j] = X.rowmajor[i][j] + Y.rowmajor[i][j];
		}
	}
	Z.copyRowToColumn();
	return Z;
}

DenseMatrix operator -(const DenseMatrix& X, const DenseMatrix& Y) {
	if (X.rows != Y.rows || X.cols != Y.cols)
	{
		//Error: mismatched dimensions in matrix substraction. Returning X
		return X;
	}
	DenseMatrix Z = DenseMatrix(X.rows,X.cols,true);
	for(int i=0; i < Z.rows; ++i)
	{
		for(int j=0; j < Z.cols; ++j)
		{
			Z.rowmajor[i][j] = X.rowmajor[i][j] - Y.rowmajor[i][j];
		}
	}
	Z.copyRowToColumn();
	return Z;
}

DenseMatrix operator *(const DenseMatrix& X, const DenseMatrix& Y) {
	if (X.cols != Y.rows)
	{
		//Error: mismatched dimensions in matrix multiplication. Returning X.
		return X;
	}
	DenseMatrix Z = DenseMatrix(X.rows,Y.cols,true);
	for(int i = 0; i < X.rows; ++i)
    	for(int j = 0; j < Y.cols; ++j)
        	for(int k = 0; k < X.cols; ++k)
            {
            	Z.rowmajor[i][j] += X.rowmajor[i][k] * Y.rowmajor[k][j];
            }
	Z.copyRowToColumn();
	return Z;
}

std::ostream& operator<<(std::ostream& os, const DenseMatrix& X)  
{  
    for (int i = 0; i < X.rows; ++i)
	{
		for (int j = 0; j < X.cols; ++j)
		{
			os << X.rowmajor[i][j] << " ";
		}
		os << endl;
	}
	return os;  
}

typedef struct LowerTriangularMatrix
{
	dtype* rowmajor;//serialized
	int dim,totalsize;
	
	LowerTriangularMatrix()
	{
		dim = 1;
		totalsize = ((dim*(dim+1))/2);
		rowmajor = new dtype[totalsize]();
	}
	
	LowerTriangularMatrix(int dim_) : dim(dim_),totalsize((dim_*(dim_+1))/2)
	{
		rowmajor = new dtype[totalsize]();
	}
	~LowerTriangularMatrix()
	{
		delete[] rowmajor;
	}
} LowerTriangularMatrix;

std::ostream& operator<<(std::ostream& os, const LowerTriangularMatrix& X)  
{  
    int ind = 0;
	for (int i = 0; i < X.dim; ++i)
	{
		for (int j = 0; j <=i; ++j)
		{
			os << X.rowmajor[ind] << " ";
			++ind;
		}
		os << endl;
	}
	return os;  
}

typedef struct Mask
{
	bool* value;
	int dim;
	std::string tmp; 
	std::vector<char> vs;
	Mask() : value(NULL),dim(0) {}
	Mask(int dim_,bool init) : dim(dim_)
	{
		value = new bool[dim]();
		for(int i=0;i<dim;++i) value[i] = init;
		tmp.reserve(dim);
		vs.reserve(dim);
	}
	~Mask()
	{
		delete[] value;
	}
} Mask;

typedef struct NNLS_Single_Input
{//Cx = b
	DenseMatrix *C;
	dtype *x,*b;
	int max_iter;
	NNLS_Single_Input(DenseMatrix* C_,dtype* b_,int max_iter_) : C(C_),b(b_),max_iter(max_iter_)
	{
		x = new dtype[C->cols]();
	}
	~NNLS_Single_Input()
	{
		delete[] x;
	}
} NNLS_Single_Input;

typedef struct NNLS_Single_State
{//Cx = b, min norm Cx-b
	DenseMatrix *C_xmask,*C_ymask;
	Mask *xmask,*infeasiblemask;
	dtype *y_masked,*x_masked,*y_masked_intermediate;
	int full_exchange_buffer,infeasible,lowest_infeasible;
	int iterations;
	bool full_exchange_mode;
	NNLS_Single_State(int rows,int cols)
	{
		xmask = new Mask(cols,false);
		infeasiblemask = new Mask(cols,false);
		C_xmask = new DenseMatrix(rows,cols,true);
		C_xmask->cols = 0;
		C_ymask = new DenseMatrix(rows,cols,true);
		x_masked = new dtype[cols]();
		y_masked = new dtype[cols]();
		y_masked_intermediate = new dtype[rows]();
		full_exchange_buffer = FULL_EXCHANGE_BUFFER;
		full_exchange_mode = true;
		lowest_infeasible = cols+1;
		iterations = 0;
	}
	~NNLS_Single_State()
	{
		delete C_xmask;
		delete C_ymask;
		delete xmask;
		delete infeasiblemask;
		delete[] x_masked;
		delete[] y_masked;
		delete[] y_masked_intermediate;
	}
} NNLS_Single_State;


#if __cplusplus > 199711L
	// typedef std::unordered_map<std::string,LowerTriangularMatrix*> CholeskyMap;
	typedef std::map<std::string,LowerTriangularMatrix*> CholeskyMap;
	// typedef std::unordered_map<std::vector<int>,LowerTriangularMatrix*> CholeskyMap2;
#else
	typedef std::map<std::string,LowerTriangularMatrix*> CholeskyMap;
	// typedef std::map<std::vector<int>,LowerTriangularMatrix*> CholeskyMap2;
#endif

inline void destroyCholeskyMap(CholeskyMap& cm)
{
	for(CholeskyMap::iterator it = cm.begin(); it != cm.end(); ++it)
	{

		delete it->second;
	}
	cm.clear();
}
// void destroyCholeskyMap2(CholeskyMap2& cm)
// {
// 	for(CholeskyMap2::iterator it = cm.begin(); it != cm.end(); ++it)
// 	{
// 		delete it->second;
// 	}
// 	cm.clear();
// }



typedef struct NNLS_Multiple_State
{//CX = B, min norm CX-B
	Mask **xmasks,**infeasiblemasks;
	int *full_exchange_buffer;
	int *lowest_infeasible;
	bool *full_exchange_mode;
	int *infeasible;
	int iterations;
	const int cols,cols_rhs;
	LowerTriangularMatrix **G;
	CholeskyMap choleskyMap;
	// CholeskyMap2 choleskyMap2;

	LowerTriangularMatrix **CFTCF;
	DenseMatrix **CGTCF;
	dtype **y_masked,**x_masked,**CGTb;

	int totalfeasible;

	NNLS_Multiple_State(int cols_, int cols_rhs_) : cols(cols_),cols_rhs(cols_rhs_)
	{
		xmasks = new Mask*[cols_rhs];
		for(int i=0;i<cols_rhs;++i) xmasks[i] = new Mask(cols,false);
		infeasiblemasks = new Mask*[cols_rhs];
		for(int i=0;i<cols_rhs;++i) infeasiblemasks[i] = new Mask(cols,false);
		
		full_exchange_buffer = new int[cols_rhs];
		for(int i=0;i<cols_rhs;++i) full_exchange_buffer[i] = FULL_EXCHANGE_BUFFER;
		lowest_infeasible = new int[cols_rhs];
		for(int i=0;i<cols_rhs;++i) lowest_infeasible[i] = cols + 1;
		full_exchange_mode = new bool[cols_rhs];
		for(int i=0;i<cols_rhs;++i) full_exchange_mode[i] = true;
		infeasible = new int[cols_rhs];
		for(int i=0;i<cols_rhs;++i) infeasible[i] = cols+1;

		CFTCF = new LowerTriangularMatrix*[cols_rhs];
		for(int i=0;i<cols_rhs;++i) CFTCF[i] = new LowerTriangularMatrix(cols);
		CGTCF = new DenseMatrix*[cols_rhs];
		for(int i=0;i<cols_rhs;++i) CGTCF[i] = new DenseMatrix(cols,cols);
		y_masked = new dtype*[cols_rhs];
		for(int i=0;i<cols_rhs;++i) y_masked[i] = new dtype[cols];
		x_masked = new dtype*[cols_rhs];
		for(int i=0;i<cols_rhs;++i) x_masked[i] = new dtype[cols];
		CGTb = new dtype*[cols_rhs];
		for(int i=0;i<cols_rhs;++i) CGTb[i] = new dtype[cols];

		G = new LowerTriangularMatrix*[cols_rhs];

		totalfeasible = 0;
	}
	~NNLS_Multiple_State()
	{
		for(int i=0;i<cols_rhs;++i)
		{
			delete CFTCF[i];
			delete CGTCF[i];
			delete xmasks[i];
			delete infeasiblemasks[i];
			delete[] x_masked[i];
			delete[] y_masked[i];
			delete[] CGTb[i];
		}

		delete[] xmasks;
		delete[] infeasiblemasks;

		delete[] full_exchange_buffer;
		delete[] lowest_infeasible;
		delete[] full_exchange_mode;
		delete[] infeasible;
		
		delete[] CFTCF;
		delete[] CGTCF;
		delete[] x_masked;
		delete[] y_masked;
		delete[] CGTb;

		delete[] G;

		destroyCholeskyMap(choleskyMap);
		// destroyCholeskyMap2(choleskyMap2);

	}
} NNLS_Multiple_State;




typedef struct NNLS_Multiple_Input
{//CX = B
	LowerTriangularMatrix *CTC;
	dtype **X,**CTB;
	const int cols_rhs;
	const int max_iter;
	bool allocateX;
	double allocate_time,init_time,
		advance_time,switch_time,determine_time,
		apply_time,normal_time,generateCGTCF_time,
		matvec_time,generateCGTb_time,vector_time,mark_time,overwrite_time;
	double choleskyCFTCF_time,choleskyFactorization_time,choleskyMapInsert_time,choleskyCheck_time,
		choleskyAllocate_time,choleskyKeyFind_time,choleskyMaskToString_time;
	NNLS_Multiple_State* state;
	NNLS_Multiple_Input(LowerTriangularMatrix* CTC_,dtype**X_,dtype**CTB_,int cols_rhs_,int max_iter_) : 
		CTC(CTC_),X(X_),CTB(CTB_),cols_rhs(cols_rhs_),max_iter(max_iter_),
		allocate_time(0),init_time(0),
		advance_time(0),switch_time(0),determine_time(0),
		apply_time(0),normal_time(0),generateCGTCF_time(0),
		matvec_time(0),generateCGTb_time(0),vector_time(0),mark_time(0),overwrite_time(0),
		choleskyCFTCF_time(0),choleskyFactorization_time(0),choleskyMapInsert_time(0),choleskyCheck_time(0),
		choleskyAllocate_time(0),choleskyKeyFind_time(0),choleskyMaskToString_time(0)
	{
		allocateX = false;
		state = new NNLS_Multiple_State(CTC->dim,cols_rhs);
	}
	NNLS_Multiple_Input(LowerTriangularMatrix* CTC_,dtype**CTB_,int cols_rhs_,int max_iter_) : 
		CTC(CTC_),CTB(CTB_),cols_rhs(cols_rhs_),max_iter(max_iter_),
		allocate_time(0),init_time(0),
		advance_time(0),switch_time(0),determine_time(0),
		apply_time(0),normal_time(0),generateCGTCF_time(0),
		matvec_time(0),generateCGTb_time(0),vector_time(0),mark_time(0),overwrite_time(0),
		choleskyCFTCF_time(0),choleskyFactorization_time(0),choleskyMapInsert_time(0),choleskyCheck_time(0),
		choleskyAllocate_time(0),choleskyKeyFind_time(0),choleskyMaskToString_time(0)
	{
		allocateX = true;
		X = new dtype*[cols_rhs];
		for(int i=0;i<cols_rhs;++i) X[i] = new dtype[CTC->dim];
		state = new NNLS_Multiple_State(CTC->dim,cols_rhs);
	}
	~NNLS_Multiple_Input()
	{
		if(allocateX)
		{
			for(int i=0;i<cols_rhs;++i) delete[] X[i];
			delete[] X;
		}
		delete state;
	}
} NNLS_Multiple_Input;







typedef struct NMF_Input
{//A: m * n, W : m * k, H: k * n
	DenseMatrix *W,*H,*A;
	int m,k,n;//data size, dimensions, data points
	int max_iter_nmf,max_iter_nnls;
	NMF_Input(DenseMatrix*W_,DenseMatrix*H_,DenseMatrix*A_,int max_iter_nmf_,int max_iter_nnls_) : 
		W(W_),H(H_),A(A_),max_iter_nmf(max_iter_nmf_),max_iter_nnls(max_iter_nnls_)
	{
		m = A->rows;
		k = W->cols;
		n = H->cols;
	}
	~NMF_Input()
	{
	}
} NMF_Input;

typedef struct NMF_State
{//A: m * n, W : m * k, H: k * n
	LowerTriangularMatrix *HHT,*WTW;
	dtype **WTA,**HAT;
	int m,k,n;//data size, dimensions, data points
	NMF_State(int m_,int k_,int n_) : m(m_),k(k_),n(n_) 
	{
		HHT = new LowerTriangularMatrix(k);
		WTW = new LowerTriangularMatrix(k);
		WTA = new dtype*[n];
		for(int col=0;col<n;++col) WTA[col] = new dtype[k];
		HAT = new dtype*[m];
		for(int col=0;col<m;++col) HAT[col] = new dtype[k];
	}
	~NMF_State()
	{
		delete HHT;
		delete WTW;
		for(int i=0;i<n;++i) delete[] WTA[i];
		delete[] WTA;
		for(int i=0;i<m;++i) delete[] HAT[i];
		delete[] HAT;
	}
} NMF_State;
