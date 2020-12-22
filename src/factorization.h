#pragma once

#include "parameters.h"
#include "vector.h"
#include "entities.h"

void cholesky_lowertriangular_cpu(LowerTriangularMatrix& G, LowerTriangularMatrix& C)
{
	int row = 0;
	int startindex = 0;
	while(startindex < C.totalsize)
	{
		int diagonalindex = startindex+row;
		dtype sqsum = 0;
		for(int col=0;col<row;++col)
		{
			int crossindex = (col*(col+1))/2;
			dtype crosssum = productsum(G.rowmajor+crossindex,G.rowmajor+startindex,col);
			int currentcolindex = startindex+col;
			int crossdiagonalindex = crossindex+col;
			G.rowmajor[currentcolindex] = 
				1.0/G.rowmajor[crossdiagonalindex]
				*(C.rowmajor[currentcolindex] - crosssum);
			sqsum += G.rowmajor[currentcolindex]*G.rowmajor[currentcolindex];
		}
		G.rowmajor[diagonalindex]=sqrt(C.rowmajor[diagonalindex]-sqsum);

		++row;
		startindex = (row*(row+1))/2;
	}
}