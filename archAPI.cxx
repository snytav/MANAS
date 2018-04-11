/*
 * archAPI.cxx
 *
 *  Created on: Apr 10, 2018
 *      Author: snytav
 */
#include<stdlib.h>
#include<string.h>

#ifndef __CUDAC__
void BlockThreadSynchronize(){}
#endif

#ifndef __CUDAC__
double MultiThreadAdd(double *address, double val)
{
    double assumed,old=*address;

#ifdef OMP_THREADS
#pragma omp critical
    {
#endif
    *address += val;

    old = *address;
#ifdef OMP_THREADS
    }
#endif

    return old;
}
#endif


const char *getErrorString(int err)
{

	return "";
}

int SetDevice(int n){return 0;}

#ifndef __CUDACC__
void AsyncCopy(double *dst,double *src,int n,int size)
{
	int j;

	for(j = 0;j < size;j++)
	{
	   dst[j] = src[j];
	}

}


int MemoryCopy(void* dst,void *src,size_t size,int dir)
{
	memcpy(dst,src,size);

	return 0;
}
#endif

#ifndef __CUDACC__
int MemoryAllocate(void** dst,size_t size)
{
	*dst = malloc(size);
    return 0;
}
#endif

#ifndef __CUDACC__
int GetDeviceMemory(size_t *m_free,size_t *m_total)
{
	*m_free = 0;
	*m_total = 0;
	return 0;
}
#endif

#ifndef __CUDACC__
int MemorySet(void *s, int c, size_t n)
{
	memset(s,c,n);
    return 0;
}
#endif


#ifndef __CUDACC__
int DeviceSynchronize()
{
    return 0;
}

 int ThreadSynchronize()
{
	 return 0;
}

 int getLastError()
{
	return 0;
}
#endif




