
void randomInit(dtype*A,int totalsize)
{
	for (int i=0;i<totalsize;++i)
	{
		A[i] = rand();
	}
}

void matvecmult_cpu(dtype*A,int rows,int cols,dtype*x,dtype*b,int factor=1)
{//A: rows * cols, x: cols, b: rows
	for(int row=0;row<rows;++row)
	{
		//interpret as row-by-column multiplication
		b[row] = 0;
		for(int col=0;col<cols;++col)
		{
			b[row] += A[row*cols + col]*x[col];
		}
		b[row] *= factor;
	}
}
void matmult_ata_cpu(dtype*C,dtype*A,int rows,int cols)
{//computes C=ATA from A.
	for(int Crow=0;Crow<cols;++Crow)
	{
		for(int Ccol=0;Ccol<cols;++Ccol)
		{
			int Cindex = Crow*cols + Ccol;
			C[Cindex] = 0;
			for(int row=0;row<rows;++row)
			{
				C[Cindex] += A[Crow + row*cols]*A[Ccol+rows*cols];
			}
		}
	}
}