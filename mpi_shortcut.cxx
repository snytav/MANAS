/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int InitMPI(int argc,char *argv[])
{
	MPI_Init(&argc,&argv);
}

int sumMPI(int size,double *jx,double *jy,double *jz)
{
    double *snd,*rcv;
    int i;
    
    snd = (double *)malloc(3*size*sizeof(double));
    rcv = (double *)malloc(3*size*sizeof(double));
    
    for(i = 0;i < size;i++)
    {   
        snd[i]          = jx[i];
        snd[i +   size] = jy[i];
        snd[i + 2*size] = jz[i];
    }

    MPI_Allreduce(snd,rcv,3*size,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
    
    for(i = 0;i < size;i++)
    {   
        jx[i] = rcv[i];
        jy[i] = rcv[i +   size];
        jz[i] = rcv[i + 2*size];
    }

    return 0;
}

int sumMPIenergy(double *e)
{
    double snd,rcv;
    
    snd = *e;

    MPI_Allreduce(&snd,&rcv,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
    
    *e = rcv;
    
    return 0;
}

int sumMPIint(int *e)
{
    int snd,rcv;

    snd = *e;

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(&snd,&rcv,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    *e = rcv;

    return 0;
}

int getRank()
{
    int rank;
    
    
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    return rank;
}

int getSize()
{
    int rank;


    MPI_Comm_size(MPI_COMM_WORLD,&rank);

    return rank;
}

int  CloseMPI()
{
	MPI_Finalize();
}

int SendComputationName(char *name)
{
	char buf[200];

	if(getSize() == 1) return 0;
//    printf("rank %d name %s \n",getRank(),name);
	strcpy(buf,name);
//    printf("rank %d buf %s \n",getRank(),buf);


    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(buf,150,MPI_CHAR,0,MPI_COMM_WORLD);
//	printf("rank %d afterbuf %s \n",getRank(),buf);
	strcpy(name,buf);

//    printf("AFTER rank %d name %s \n",getRank(),name);

    return 0;
}
