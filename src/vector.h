#pragma once

#include "parameters.h"

dtype vectordot(dtype*a,dtype*b,int size)
{
	dtype sum = 0;
	for(int i=0;i<size;++i) sum += a[i]*b[i];
	return sum;
}
dtype productsum(dtype*a,dtype*b,int size)
{
	dtype sum=0;
	for(int i=0;i<size;++i) sum += a[i]*b[i];
	return sum;
}
void vectoradd(dtype*a,dtype*b,int size,dtype factor=1)
{
	for(int i=0;i<size;++i)
	{
		b[i] += factor*a[i];
	}
}
void vectorsub(dtype*a,dtype*b,int size)
{
	for(int i=0;i<size;++i) b[i] -= a[i];
}
void vectornegate(dtype*a,dtype*b,int size)
{
	for(int i=0;i<size;++i) b[i] = -a[i];
}
void vectorinit(dtype*a,int size,dtype value)
{
	for(int i=0;i<size;++i) a[i] = value;
}