#pragma once

#include <sstream>

#include "entities.h"

void applyColumnMask(DenseMatrix& original,DenseMatrix& masked,Mask& mask,bool toggle=false)
{
	masked.cols = 0;
	masked.totalsize = 0;
	for(int i=0;i<mask.dim;++i)
	{
		if( mask.value[i] xor toggle )
		{
			// printf("masked index: %d\n",i);
			masked.colmajor[masked.cols] = original.colmajor[i];
			++masked.cols;
		}
	}
	masked.totalsize = masked.cols*masked.rows;
}
void applyVectorMask(dtype*original,dtype*masked,Mask& mask,bool toggle=false)
{
	int counter = 0;
	for(int i=0;i<mask.dim;++i)
	{
		if( mask.value[i] xor toggle )
		{
			masked[counter] = original[i];
			++counter;
		}
	}
}
void overwriteOriginalWithMask(dtype*original,dtype*masked,Mask& mask)
{
	int counter = 0;
	for(int i=0;i<mask.dim;++i)
	{
		if( mask.value[i] )
		{
			original[i] = masked[counter];
			++counter;
		}
		else
		{
			original[i] = 0;
		} 
	}
}
int markInfeasible(Mask& infeasiblemask,dtype*x_masked,dtype*y_masked,Mask& xmask)
{
	int xcounter = 0;
	int ycounter = 0;
	int infeasible = 0;
	for(int i=0;i<xmask.dim;++i)
	{
		if( xmask.value[i] )
		{
			infeasiblemask.value[i] = x_masked[xcounter] < 0;
			++xcounter;
		}
		else
		{
			infeasiblemask.value[i] = y_masked[ycounter] < 0;
			++ycounter;
		}

		if ( infeasiblemask.value[i] ) { ++infeasible; }
	}
	return infeasible;
}

void switchSets(Mask& infeasiblemask,Mask& xmask,bool bpp)
{
	for(int i=xmask.dim-1;i>=0;--i)
	{
		if( infeasiblemask.value[i] )
		{
			xmask.value[i] = not xmask.value[i];
			if( not bpp ) return;
		}
	}
}

void generateCGTCF(DenseMatrix& CGTCF,LowerTriangularMatrix& CTC,Mask& xmask)
{
	int rowmap[xmask.dim];
	int colmap[xmask.dim];
	int newrows = 0;
	for(int row = 0; row < xmask.dim; ++row)
	{
		if( not xmask.value[row] )
		{
			rowmap[newrows] = row;
			++newrows;
		}
	}
	CGTCF.rows = newrows;
	int newcols = 0;
	for(int col = 0; col < xmask.dim; ++col)
	{
		if( xmask.value[col] )
		{
			colmap[newcols] = col;
			++newcols;
		}
	}
	CGTCF.cols = newcols;

	CGTCF.totalsize = CGTCF.rows*CGTCF.cols;

	for(int row = 0; row < CGTCF.rows; ++row)
	{
		for(int col = 0; col < CGTCF.cols; ++col )
		{
			int oldrow = rowmap[row];
			int oldcol = colmap[col];
			int index = (oldcol > oldrow) ? ((oldcol*(oldcol+1))/2+oldrow) : ((oldrow*(oldrow+1))/2+oldcol);
			CGTCF.colmajor[col][row] = CTC.rowmajor[index];
		}
	}
}
void generateCFTCF(LowerTriangularMatrix& CFTCF,LowerTriangularMatrix& CTC,Mask& columnmask)
{
	int colmap[columnmask.dim];
	int newdim = 0;
	for(int i = 0; i < columnmask.dim; ++i)
	{
		if( columnmask.value[i] )
		{
			colmap[newdim] = i;
			++newdim;
		}
	}
	CFTCF.dim = newdim;
	CFTCF.totalsize = (newdim*(newdim+1))/2;

	for(int row = 0; row < CFTCF.dim; ++row)
	{
		int index = (row*(row+1))/2;
		for(int col = 0; col <= row; ++col )
		{
			int oldrow = colmap[row];
			int oldcol = colmap[col];
			int oldindex = (oldcol > oldrow) ? ((oldcol*(oldcol+1))/2+oldrow) : ((oldrow*(oldrow+1))/2+oldcol);
			// std::cout << "oldcol,oldrow,oldindex: " << oldcol<<","<<oldrow<<","<<oldindex<<std::endl;
			// std::cout << "Press once" << std::endl;
			// std::cin.ignore();
			//std::cout << "CTC totalsize: " << CTC.totalsize << std::endl;
			// std::cout << "pointer: " << CTC.rowmajor << std::endl;
			CFTCF.rowmajor[index + col] = CTC.rowmajor[oldindex];
			// std::cout << "Press twice" << std::endl;
			// std::cin.ignore();
		}
	}
}
void generateCGTb(dtype* CTb,dtype* CGTb,Mask& xmask)
{
	int newdim = 0;
	for(int i = 0; i < xmask.dim; ++i)
	{
		if( not xmask.value[i] )
		{
			CGTb[newdim] = CTb[i];
			++newdim;
		}
	}
}
std::string maskToString(Mask& mask)
{
	std::stringstream stream;
	// std::string key = "";
	for(int i=0;i<mask.dim;++i)
	{
		if( mask.value[i] )
		{
			stream << i << "_";
		}
	}
	return stream.str();
}
void maskToString3(Mask& mask)
{
	std::stringstream stream;
	mask.tmp = "";
	for(int i=0;i<mask.dim;++i)
	{
		if( mask.value[i] )
		{
			stream << i << "_";
		}
	}
	mask.tmp = stream.str();
}
void num2string_reverse(int num,std::string& str)
{
	while(num > 0)
	{
		str += (char)((num%10) + 48);
		num/=10;
	}
}
void maskToString4(Mask& mask)
{
	mask.tmp = "";
	for(int i=0;i<mask.dim;++i)
	{
		if( mask.value[i] )
		{
			num2string_reverse(i,mask.tmp);
			mask.tmp+='_';
		}
	}
}
void num2vec_reverse(int num,std::vector<char>& vec)
{
	while(num > 0)
	{
		vec.push_back( (char)((num%10) + 48) );
		num/=10;
	}
}
void maskToString5(Mask& mask)
{
	mask.vs.clear();
	for(int i=0;i<mask.dim;++i)
	{
		if( mask.value[i] )
		{
			num2vec_reverse(i,mask.vs);
			mask.vs.push_back('_');
		}
	}
	mask.vs.push_back('\0');
}
// void maskToString2(Mask& mask,std::string& key)
// {
// 	// std::stringstream stream;
// 	key.reserve(5*mask.dim);
// 	key = "";
// 	for(int i=0;i<mask.dim;++i)
// 	{
// 		if( mask.value[i] )
// 		{
// 			// stream << i;
// 			// key += stream.str() + "_";
// 			// stream.str(std::string());
// 			char num[11];
// 			key += itoa(i,num,16);
// 			key += "_";
// 			// sprintf(num,"%d",i);
// 			// key += num;
// 			// key.append(itoa(i,num,36));
// 			// key.append("_");
// 		}
// 	}
// 	//key = key.substr(0, key.size()-1);
// }

void maskToVector(Mask& mask,std::vector<int>& key)
{
	for(int i=0;i<mask.dim;++i)
	{
		if( mask.value[i] )
		{
			key.push_back(i);
		}
	}
}

