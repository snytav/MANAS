/*
 * gpu_plasma.h
 *
 *  Created on: Aug 21, 2013
 *      Author: snytav
 */
//#include "cuPrintf.cu"

// Plans 10th November:


//        3.4. memory copies in GPUCell
// 4. Derivative class for Xeon
// 5. Derivative class for OpenMP
//  To be completed to Nov 14th !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// 6. performance comparison
// 7. Launch at Phi (accuracy check)
// 8. Many-node performance tests with Phi



#ifndef GPU_PLASMA_H_
#define GPU_PLASMA_H_
#include<stdlib.h>
#include<stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>

#include <sys/types.h>
#include <sys/stat.h>

//#include <unistd.h>
//#include <stdio.h>
#include <errno.h>

#ifdef __CUDACC__
#include <nvToolsExtCuda.h>
#include <nvToolsExtCudaRt.h>
#endif


#include <time.h>

//#ifdef __OMP__
#include <omp.h>
//#endif

#ifdef __CUDACC__
#include <cuda.h>
#endif

#include "archAPI.h"
#include "rnd.h"
#include "plasma.h"
#include "gpucell.h"
#include "mpi_shortcut.h"

#include <sys/resource.h>
#include <stdint.h>

#include <sys/sysinfo.h>
#include <sys/time.h>

//#include "read_particles.c"

#include "init.h"
#include "diagnose.h"

#include<string>
#include <iostream>

#include "particle_target.h"

#include "params.h"

#include "memory_control.h"



using namespace std;


double get_meminfo(void)
{
	FILE *f;
	char str[100];
	int  mem_free;
	double dmem;
   // return 0.0;

	system("free>&free_mem_out.dat");


	if((f = fopen("free_mem_out.dat","rt")) == NULL) return 0.0;

	fgets(str,100,f);
	fgets(str,100,f);

	mem_free = atoi(str + 30);

	dmem = (((double)mem_free)/1024)/1024;

	return dmem;

}

double get_meminfo1(void)
{
	double retval=0;
	char tmp[256]={0x0};
	/* note= add a path to meminfo like /usr/bin/meminfo
	   to match where meminfo lives on your system */
	FILE *shellcommand=popen("meminfo","r");
	while(fgets(tmp,sizeof(tmp),shellcommand)!=NULL)
	{
		if(memcmp(tmp,"Mem:",4)==0)
		{
			int	wordcount=0;
			char delimiter[10];

			strcpy(delimiter," ");

			char *p=strtok(tmp,delimiter);
			while(*p)
			{
				wordcount++;
				if(wordcount==3) retval=atof(p);
			}
		}
	}
	pclose(shellcommand);
	return retval;
}




#define FORTRAN_ORDER




//#else
//double cuda_atomicAdd(double *address, double val)
//{
//    double assumed,old=*address;
//
//    old += val;
//
//    return old;
//}
//#endif


const int flagCPUandGPUrun = 1;



//template <template <class Particle> class Cell >
//#ifdef __CUDACC__
//__global__
//#endif
//void printParticle(Cell<Particle>  **cells,int num,int sort)
//{
//	unsigned int nx = blockIdx.x;
//	unsigned int ny = blockIdx.y;
//	unsigned int nz = blockIdx.z;
//
//	Cell<Particle>  *c,*c0 = cells[0],nc;
//	__shared__ extern CellDouble fd[9];
//	Particle p;
//
//	c = cells[ c0->getGlobalCellNumber(nx,ny,nz)];
//
//	nc = *c;
//    if(nc.number_of_particles < threadIdx.x) return;
//
//	nc.readParticleFromSurfaceDevice(threadIdx.x,&p);
//
//		if(p.fortran_number == num && (int)p.sort == sort)
//		{
//			printf("particle-print %5d thread %3d cell (%d,%d,%d) sort %d  %25.15e,%25.15e,%25.15e \n",p.fortran_number,threadIdx.x,c->i,c->l,c->k,(int)p.sort,p.x,p.y,p.z);
//		}
//}


//template<int dims>
typedef void (*SingleNodeFunctionType)(GPUCell<Particle,DIMENSIONS>  **cells,KernelParams *params,
		unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		unsigned int nx,unsigned int ny,unsigned int nz
		);
//typedef void (*SingleNodeFunctionType)(void  *,unsigned int,unsigned int,unsigned int);







//template <template <class Particle,int dims> class Cell,int dims >
//__device__
//void emh2(
//		 Cell<Particle,dims>  **cells,
//		 KernelParams *params,
//   		 unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
//		 unsigned int tnx,unsigned int tny,unsigned int tnz
//		)
//{
//	unsigned int nx = bk_nx;
//	unsigned int ny = bk_ny;
//	unsigned int nz = bk_nz;
//	Cell<Particle,dims>  *c0 = cells[0];
//
//	emh2_Element(c0,params->i_s+nx,params->l_s+ny,params->k_s+nz,params->Q,params->H);
//}

template <template <class Particle,int dims> class Cell,int dims >
__device__
void emh2(
		 Cell<Particle,dims>  **cells,
		 KernelParams *params,
   		 unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		 unsigned int tnx,unsigned int tny,unsigned int tnz
		)
{
	unsigned int nx = bk_nx;
	unsigned int ny = bk_ny;
	unsigned int nz = bk_nz;
	Cell<Particle,dims>  *c0 = cells[0];

	emh2_Element(c0,params->i_s+nx,params->l_s+ny,params->k_s+nz,params->Q,params->H);
}

template <template <class Particle,int dims> class Cell,int dims>
__device__
void GPU_add(
		 Cell<Particle,dims>  **cells,
				 KernelParams *params,
		   		 unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
				 unsigned int tnx,unsigned int tny,unsigned int tnz
		)
{
	unsigned int nx = bk_nx;
		unsigned int ny = bk_ny;
		unsigned int nz = bk_nz;
		Cell<Particle,dims>  *c0 = cells[0];

//	printf("bk_nx %d bk_ny %d bk_nz %d tnx %d tny %d tnz %d \n",bk_nx,bk_ny,bk_nz,tnx,tny,tnz);
	add_Element(c0,params->i_s+nx,params->l_s+ny,params->k_s+nz,params->B0,params->H);

}

template <template <class Particle,int dims> class Cell,int dims >
__device__
void MagneticField_SingleNode(
		 Cell<Particle,dims>  **cells,
		 KernelParams *params,
		 unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		 unsigned int tnx,unsigned int tny,unsigned int tnz
		)
{
	unsigned int nx = bk_nx;
	unsigned int ny = bk_ny;
	unsigned int nz = bk_nz;
	Cell<Particle,dims>  *c0 = cells[0];

	emh1_Element(c0,params->i_s+nx,params->l_s+ny,params->k_s+nz,
			params->Q,params->H,params->E1,params->E2,
			params->c1,params->c2,params->dx1,params->dy1,
			params->dz1,params->dx2,params->dy2,params->dz2);
}




template <template <class Particle,int dims> class Cell,int dims >
__device__ void GPU_SetAllCurrentsToZero_SingleNode(Cell<Particle,dims>  **cells,KernelParams *params,
		unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		unsigned int nx,unsigned int ny,unsigned int nz)
{
	Cell<Particle,dims>  *c,*c0 = cells[0],nc;
	__shared__ extern CellDouble fd[9];

//	Cell<Particle,dims>  **cells = (Cell<Particle,dims>  **)ptr;

	c = cells[ c0->getGlobalCellNumber(bk_nx,bk_ny,bk_nz)];

	nc = *c;

	nc.SetAllCurrentsToZero(nx,ny,nz);
}

__device__ void writeCurrentComponent(CellDouble *J,CurrentTensorComponent *t1,CurrentTensorComponent *t2,int pqr2)
{
    MultiThreadAdd (J->getp(t1->i11,t1->i12,t1->i13),t1->t[0]);
    MultiThreadAdd (J->getp(t1->i21,t1->i22,t1->i23),t1->t[1]);
    MultiThreadAdd (J->getp(t1->i31,t1->i32,t1->i33),t1->t[2]);
    MultiThreadAdd (J->getp(t1->i41,t1->i42,t1->i43),t1->t[3]);

    if(pqr2 == 2)
    {
        MultiThreadAdd(J->getp(t2->i11,t2->i12,t2->i13),t2->t[0]);
        MultiThreadAdd(J->getp(t2->i21,t2->i22,t2->i23),t2->t[1]);
        MultiThreadAdd(J->getp(t2->i31,t2->i32,t2->i33),t2->t[2]);
        MultiThreadAdd(J->getp(t2->i41,t2->i42,t2->i43),t2->t[3]);
    }

}

template <template <class Particle,int dims> class Cell,int dims >
__device__ void CurrentPeriodic(Cell<Particle,dims>  **cells,
		                        KernelParams *params,
						        unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
						        unsigned int tnx,unsigned int tny,unsigned int tnz
                             )
{
	unsigned int nx = bk_nx;
	unsigned int nz = bk_nz;
	Cell<Particle,dims>  *c0 = cells[0];

	periodicCurrentElement(c0,nx+params->i_s,nz+params->k_s,params->E,
			               params->dir,params->dirE,params->N);
}

template <template <class Particle,int dims> class Cell,int dims >
__device__ void GPU_periodic_SingleNode(Cell<Particle,dims>  **cells,
		                    KernelParams *params,
				            unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
				            unsigned int tnx,unsigned int tny,unsigned int tnz
                             )
{
	unsigned int nx = bk_nx;
	//unsigned int ny = blockIdx.y;
	unsigned int nz = bk_nz;
	Cell<Particle,dims>  *c0 = cells[0];

	periodicElement(c0,nx+params->i_s,nz+params->k_s,params->E, params->dir,
			params->to,params->from);

}

template <template <class Particle,int dims> class Cell,int dims>
__device__ void GPU_eme_SingleNode(

		            Cell<Particle,dims>  **cells,
		            KernelParams *params,
		            unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		            unsigned int tnx,unsigned int tny,unsigned int tnz
		)
{
	unsigned int nx = bk_nx*params->blockDim_x + tnx;
	unsigned int ny = bk_ny*params->blockDim_y + tny;
	unsigned int nz = bk_nz*params->blockDim_z + tnz;
	Cell<Particle,dims>  *c0 = cells[0];

	emeElement(c0,params->i_s+nx,params->l_s+ny,params->k_s+nz,params->E,params->H1,params->H2,
			params->J,params->c1,params->c2,params->tau,
			params->dx1,params->dy1,params->dz1,
			params->dx2,params->dy2,params->dz2);
}

template <template <class Particle,int dims> class Cell,int dims >
__device__ void GPU_MakeDepartureLists_SingleNode(Cell<Particle,dims>  **cells,
		                    KernelParams *params,
				            unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
				            unsigned int tnx,unsigned int tny,unsigned int tnz
		    )
{


	    unsigned int nx = bk_nx;
		unsigned int ny = bk_ny;
		unsigned int nz = bk_nz;
		int ix,iy,iz,n;

		Particle p;
		Cell<Particle,dims>  *c,*c0 = cells[0],nc,*new_c;
		c = cells[ n = c0->getGlobalCellNumber(nx,ny,nz)];

		c->departureListLength = 0;
		for(ix = 0;ix < 3;ix++)
		{
			for(iy = 0;iy < 3;iy++)
			{
				for(iz = 0;iz < 3;iz++)
				{
					c->departure[ix][iy][iz]      = 0;
				}
			}
		}
		c->departureListLength  = 0;

		for(int num = 0;num < c->number_of_particles; num++)
			{
			c->readParticleFromSurfaceDevice(num,&p);

				if(!c->isPointInCell(p.GetX()))   //check Paricle = operator !!!!!!!!!!!!!!!!!!!!!!!!!!!
				{
					c->removeParticleFromSurfaceDevice(num,&p,&(c->number_of_particles));
					c->flyDirection(&p,&ix,&iy,&iz);
//					if(p.fortran_number == 325041 && p.sort == 2) {
//						d_stage[0] = ix;
//						d_stage[1] = iy;
//						d_stage[2] = iz;
//					}


                    if(c->departureListLength == PARTICLES_FLYING_ONE_DIRECTION)
                    {
                    	params->d_stage[0] = TOO_MANY_PARTICLES;
                    	params->d_stage[1] = c->cnum.x;
                    	params->d_stage[2] = c->cnum.y;
                    	params->d_stage[3] = c->cnum.z();
                    	params->d_stage[4] = ix;
                    	params->d_stage[5] = iy;
                    	params->d_stage[6] = iz;
                    	return;
                    }
					c->departureListLength++;
					int num1 = c->departure[ix][iy][iz];

					c->departureList[ix][iy][iz][num1] = p;
//					if(p.fortran_number == 325041 && p.sort == 2) {
//						d_stage[4] = num1;
//						d_stage[5] = c->departureList[ix][iy][iz][num1].fortran_number;
//
//					}

					c->departure[ix][iy][iz] += 1;
					num--;
				}
			}
}


template <template <class Particle,int dims> class Cell,int dims>
__device__ void ArrangeFlights(Cell<Particle,dims>  **cells,
		KernelParams *params,
        unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		unsigned int tnx,unsigned int tny,unsigned int tnz
		)
{
	unsigned int nx = bk_nx;
	unsigned int ny = bk_ny;
	unsigned int nz = bk_nz;
	int ix,iy,iz,snd_ix,snd_iy,snd_iz,num,pos,n;
	Particle p;

	Cell<Particle,dims>  *c,*c0 = cells[0],nc,*snd_c;

		c = cells[ n = c0->getGlobalCellNumber(nx,ny,nz)];

		for(ix = 0;ix < 3;ix++)
			for(iy = 0;iy < 3;iy++)
				for(iz = 0;iz < 3;iz++)
				{
					int index = ix*9 +iy*3 +iz;
					n = c0->getWrapCellNumber(nx+ix-1,ny+iy-1,nz+iz-1);

		            snd_c  = cells[ n ];
					if(nx == 24 && ny == 2 && nz == 2)
					{

						params->d_stage[index*4]   = snd_c->cnum.x;
						params->d_stage[index*4+1] = snd_c->cnum.y;
						params->d_stage[index*4+2] = snd_c->cnum.z();
						params->d_stage[index*4+3] = snd_c->departureListLength;
					}

					snd_ix = ix;
					snd_iy = iy;
					snd_iz = iz;

					c->inverseDirection(&snd_ix,&snd_iy,&snd_iz);

					num = snd_c->departure[snd_ix][snd_iy][snd_iz];
					if(c->cnum.x == 1 && (snd_c->departure[snd_ix][snd_iy][snd_iz] > 0))
					{

#ifdef TRANSFER_OUTPUT
						printf("basic cell (%d,%d,%d) ARRANGE cell:(%d,%d) dir:(%d,%d,%d) %d\n",
								c->cnum.x,c->cnum.y,c->cnum.z(),
								snd_c->cnum.y,snd_c->cnum.z(),
								snd_ix,snd_iy,snd_iz,
								snd_c->departure[snd_ix][snd_iy][snd_iz]);
#endif
					}

					for(int i = 0;i < num;i++)
					{
						p = snd_c->departureList[snd_ix][snd_iy][snd_iz][i];
						c->Insert(p);

					}


				}

}


template <template <class Particle,int dims> class Cell,int dims>
__device__ void AddParticlesToCells(Cell<Particle,dims>  **cells,
		KernelParams *params,
        unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		unsigned int tnx,unsigned int tny,unsigned int tnz
		)
{
	unsigned int nx = bk_nx;
	unsigned int ny = bk_ny;
	unsigned int nz = bk_nz;//,i,l,k;

	Cell<Particle,dims>  *c0 = cells[0],nc,*c;
//	int n  = c0->getGlobalCellNumber(i,l,k);

	c = cells[ c0->getGlobalCellNumber(nx,ny,nz)];
	nc = *c;
	if(nc.cnum.x > 0    && nc.cnum.x < nc.mesh.x - 1 &&
	   nc.cnum.y > 0    && nc.cnum.y < nc.mesh.y - 1 &&
	   nc.cnum.z() > 0  && nc.cnum.z() < nc.mesh.z() - 1
	  )
	{
		return; //NOT A BOUNDARY CELL
	}

	if(nc.cnum.x == 0)
	{
		int qq = 0;
	}
    Particle *p_array = (Particle*)params->d_p_send_array;

    for(int i = 0;i < params->p_send_array_size;i++)
    {
    	Particle *p = &(p_array[i]);

    	if(c->isPointInCell(p->GetX()))
    	{
    	   if( (c->Insert(*p) == true) )
    	   {
#ifdef TRANSFER_OUTPUT
    	   	  printf("particle %d (%e,%e,%e) is number %d in cell (%d,%d,%d) ON DEVICE\n",
    	   			 i+1,
    	   			 p->X.x,p->X.y,p->X.z(),c->number_of_particles,c->cnum.x,c->cnum.y,c->cnum.z());
#endif
    	   }
    	}
    }



}


template <template <class Particle,int dims> class Cell,int dims>
__device__ void GPU_getCellEnergy_SingleNode(Cell<Particle,dims>  **cells,KernelParams *params,
		unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		unsigned int nx,unsigned int ny,unsigned int nz)
{
	unsigned int i = bk_nx;
	unsigned int l= bk_ny;
	unsigned int k = bk_nz;
	Cell<Particle,dims>  *c0 = cells[0],nc;
	double t,ex,ey,ez;
	__shared__ extern CellDouble fd[9];
	int n  = c0->getGlobalCellNumber(i,l,k);

	ex = params->d_Ex[n];
	ey = params->d_Ey[n];
	ez = params->d_Ez[n];

	t = ex*ex+ey*ey+ez*ez;

	MultiThreadAdd(&(params->d_ee),t);
}

template <template <class Particle,int dims> class Cell,int dims>
__device__ void GPU_SetFieldsToCells_SingleNode(Cell<Particle,dims>  **cells,
		KernelParams *params,
		unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		unsigned int tnx,unsigned int tny,unsigned int tnz
		)
{
	unsigned int nx = bk_nx;
	unsigned int ny = bk_ny;
	unsigned int nz = bk_nz;
	Cell<Particle,dims>  *c,*c0 = cells[0],nc;

	c = cells[ c0->getGlobalCellNumber(nx,ny,nz)];

	nc = *c;

	nc.readFieldsFromArrays(params->d_Ex,params->d_Ey,params->d_Ez,
			params->d_Hx,params->d_Hy,params->d_Hz,tnx,tny,tnz);
}

template <template <class Particle,int dims> class Cell,int dims>
__device__ void GPU_WriteAllCurrents_SingleNode(Cell<Particle,dims>  **cells,
		KernelParams *params,
		unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		unsigned int tnx,unsigned int tny,unsigned int tnz
		)
{
	unsigned int nx = bk_nx;
	unsigned int ny = bk_ny;
	unsigned int nz = bk_nz,n_global;

	Cell<Particle,dims>  *c,*c0 = cells[0],nc;
	__shared__ extern CellDouble fd[9];

	c = cells[ c0->getGlobalCellNumber(nx,ny,nz)];


	 nc = *c;





	             int i1,l1,k1;
	        	 i1 = tnx;
	        	 l1 = tny;
	        	 k1 = tnz;
    	         int n = nc.getFortranCellNumber(
    	        		 nc.cnum.x  +i1-1,
    	        		 nc.cnum.y  +l1-1,
    	        		 nc.cnum.z()+k1-1);

    	         if (n < 0 ) n = -n;
        		 double t,t_x,t_y;//,jx_p,jy_p;
		         //int i_f,l_f;//k_f;
		         t_x = nc.Jx->get(i1,l1,k1);
		         int3 i3 = nc.getCellTripletNumber(n);

		         double dbg_t = params->d_Jx[n];

		         MultiThreadAdd(&(params->d_Jx[n]),t_x);
#ifdef CURRENT_TENSOR_MOVE_PRINTS
		         printf("%3d %3d %3d n %5d nt %d rank %d c (%3d,%3d,%3d) global (%5d,%5d,%5d) jxb %25.15e + %25.15e jxa %25.15e\n",
		        		  nc.cnum.x,nc.cnum.y,nc.cnum.z(),n,
		        		  params->nt,getRank(),
		        		  i1,l1,k1,
		        		  nc.cnum.x  +i1-1,nc.cnum.y  +l1-1,nc.cnum.z()+k1-1,
		        		  dbg_t,
		        		  t_x,
		        		  params->d_Jx[n]
		        	   );
#endif

		         t_y= nc.Jy->get(i1,l1,k1);
		         MultiThreadAdd(&(params->d_Jy[n]),t_y);
		         t = nc.Jz->get(i1,l1,k1);
		         MultiThreadAdd(&(params->d_Jz[n]),t);

}


template <template <class Particle,int dims> class Cell,int dims>
__device__ void GPU_GetCellNumbers_SingleNode(Cell<Particle,dims>  **cells,
		KernelParams *params,
		unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		unsigned int tnx,unsigned int tny,unsigned int tnz
		)
{
		Cell<Particle,dims>  *c;//,nc;
		c = cells[bk_nx];

		params->numbers[bk_nx] = (*c).number_of_particles;
}

template <template <class Particle,int dims> class Cell,int dims>
__device__ void GPU_GetFlownBeamNumber_SingleNode(Cell<Particle,dims>  **cells,
		KernelParams *params,
		unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		unsigned int tnx,unsigned int tny,unsigned int tnz
		)
{
	    Cell<Particle,dims>  *c,*c0 = cells[0];
	    int cell_num;
		c = cells[ cell_num =  c0->getGlobalCellNumber(bk_nx,bk_ny,bk_nz)];

		if(c->cnum.x == 19 && c->cnum.y == 18 && c->cnum.z() == 0)
		{
					int qq  = 0;
				}

		printf("(%d,%d,%d) FLOwN BEAM %e\n",c->cnum.x,c->cnum.y,c->cnum.z(),(*c).beam_boundary_counter);

		if((*c).beam_boundary_counter > 0)
		{
			int qq = 0;
		}

		MultiThreadAdd(&(params->flown_beam_particles),(*c).beam_boundary_counter);
}

template <template <class Particle,int dims> class Cell,int dims>
__device__ void GPU_StepAllCells_SingleNode(Cell<Particle,dims>  **cells,
		KernelParams *params,
		unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		unsigned int tnx,unsigned int tny,unsigned int tnz
		                         )
{

	int index  = tnx;
	Cell<Particle,dims>  *c,*c0 = cells[0];
	c = cells[ c0->getGlobalCellNumber(bk_nx,bk_ny,bk_nz)];

	c->beam_boundary_counter = 0;

#ifdef __CUDACC__
	__shared__ extern
#endif
	CellDouble fd[9];
	CellDouble *c_jx,*c_jy,*c_jz,*c_ex,*c_ey,*c_ez,*c_hx,*c_hy,*c_hz;
	CurrentTensor t1,t2;
	int pqr2;
	Particle p;



	c_ex = &(fd[0]);
	c_ey = &(fd[1]);
	c_ez = &(fd[2]);

	c_hx = &(fd[3]);
	c_hy = &(fd[4]);
	c_hz = &(fd[5]);

	c_jx = &(fd[6]);
	c_jy = &(fd[7]);
	c_jz = &(fd[8]);


#ifdef __CUDACC__
	while(index < CellExtent*CellExtent*CellExtent)
	{

		copyCellDouble(c_ex,c->Ex,index);
		copyCellDouble(c_ey,c->Ey,index);
		copyCellDouble(c_ez,c->Ez,index);

		copyCellDouble(c_hx,c->Hx,index);
		copyCellDouble(c_hy,c->Hy,index);
		copyCellDouble(c_hz,c->Hz,index);

		copyCellDouble(c_jx,c->Jx,index);
		copyCellDouble(c_jy,c->Jy,index);
		copyCellDouble(c_jz,c->Jz,index);

		index += params->particles_processed_by_a_single_thread;
	}
#else
	int index1 = 0;
	while(index1 < CellExtent*CellExtent*CellExtent)
		{

			copyCellDouble(c_ex,c->Ex,index1);
			copyCellDouble(c_ey,c->Ey,index1);
			copyCellDouble(c_ez,c->Ez,index1);

			copyCellDouble(c_hx,c->Hx,index1);
			copyCellDouble(c_hy,c->Hy,index1);
			copyCellDouble(c_hz,c->Hz,index1);

			copyCellDouble(c_jx,c->Jx,index1);
			copyCellDouble(c_jy,c->Jy,index1);
			copyCellDouble(c_jz,c->Jz,index1);

			double *d_dst,*d_src;

			d_dst = (double *)c_jx;
			d_src = (double *)(c->Jx);

            int n = index1;

			d_dst[n] = d_src[n];


			index1 ++;
		}
#endif
BlockThreadSynchronize();
    index = tnx;
	while(index < c->number_of_particles)
	{
        c->Move (index,&pqr2,&t1,&t2,
        		params->mass,
        		params->q_mass,
        		params->d_ctrlParticles,
        		params->jmp,
        		c_ex,c_ey,c_ez,c_hx,c_hy,c_hz,&(params->flown_beam_particles),params->nt
                );

        c->readParticleFromSurfaceDevice(index,&p);
#ifdef CURRENT_TENSOR_MOVE_PRINTS
        printf("Jx t1 %10d %5d %5d %3d rank %d nt %d sort %d fn %10d %d %d %d %25.15e \n",c->cnum.x,c->cnum.y,c->cnum.z(),index,getRank(),params->nt,p.sort,p.fortran_number,
                		                                               t1.Jx.i11, t1.Jx.i12, t1.Jx.i13,t1.Jx.t[0]);
        printf("Jx t2 %10d %5d %5d %3d rank %d nt %d sort %d fn %10d %d %d %d %25.15e \n",c->cnum.x,c->cnum.y,c->cnum.z(),index,getRank(),params->nt,p.sort,p.fortran_number,
                		                                               t1.Jx.i21, t1.Jx.i22, t1.Jx.i23,t1.Jx.t[1]);
        printf("Jx t3 %10d %5d %5d %3d rank %d nt %d sort %d fn %10d %d %d %d %25.15e \n",c->cnum.x,c->cnum.y,c->cnum.z(),index,getRank(),params->nt,p.sort,p.fortran_number,
                		                                               t1.Jx.i31, t1.Jx.i32, t1.Jx.i33,t1.Jx.t[2]);
        printf("Jx t2 %10d %5d %5d %3d rank %d nt %d sort %d fn %10d %d %d %d %25.15e \n",c->cnum.x,c->cnum.y,c->cnum.z(),index,getRank(),params->nt,p.sort,p.fortran_number,
                		                                               t1.Jx.i41, t1.Jx.i42, t1.Jx.i43,t1.Jx.t[3]);
#endif

        writeCurrentComponent(c_jx,&(t1.Jx),&(t2.Jx),pqr2);

        writeCurrentComponent(c_jy,&(t1.Jy),&(t2.Jy),pqr2);
        writeCurrentComponent(c_jz,&(t1.Jz),&(t2.Jz),pqr2);

        index += params->particles_processed_by_a_single_thread;
	}


    BlockThreadSynchronize();

    index= tnx;
    AsyncCopy((double*)c->Jx,(double*)c_jx,index,CellExtent*CellExtent*CellExtent);

#ifdef __CUDACC__
	while(index < CellExtent*CellExtent*CellExtent)
	{
//      	copyCellDouble(c->Jx,c_jx,index);
    	copyCellDouble(c->Jy,c_jy,index);
    	copyCellDouble(c->Jz,c_jz,index);

    	index += params->particles_processed_by_a_single_thread;
    }
#else
	index1 = 0;
	while(index1 < CellExtent*CellExtent*CellExtent)
		{
//	      	copyCellDouble(c->Jx,c_jx,index1);
	    	copyCellDouble(c->Jy,c_jy,index1);
	    	copyCellDouble(c->Jz,c_jz,index1);

	    	index1 ++;
	    }

#endif

    BlockThreadSynchronize();



    c->busyParticleArray = 0;

}

template <template <class Particle,int dims> class Cell,int dims>
#ifdef __CUDACC__
__device__
#endif
void GPU_Diagnose_SingleNode(Cell<Particle,dims>  **cells,
		KernelParams *params,
		unsigned int bk_nx,unsigned int bk_ny,unsigned int bk_nz,
		unsigned int tnx,unsigned int tny,unsigned int tnz
		                         )
{
	float *d;
    int cnp,sts;
	int index  = tnx;
	Cell<Particle,dims>  *c,*c0 = cells[0];
	Particle p;
	ParticleFloatArraysGroup *pag;
	ParticleFloatArrays      *pa;
	c = cells[ c0->getGlobalCellNumber(bk_nx,bk_ny,bk_nz)];

//	printf("%d %d %d num %d \n",c->cnum.x,c->cnum.y,c->cnum.z(), c->number_of_particles);

//	printf("GPU_Diagnose_SingleNode\n");
//	return;



	pag = params->d_part_diag;



	for(int i = 0;i < c->number_of_particles;i++)
	{
	    c->readParticleFromSurfaceDevice(i,&p);
	    int sort = p.sort;
	    pa = &((*pag)[sort]);

	    if(sort == 0 && p.fortran_number == 1)
	    {
	    	printf("QQ\n");
	    	printf("sort %d fn %10d x %15.5e \n",p.sort,p.fortran_number,p.X.x);
	    }

        pa->dbg_x[p.fortran_number-1]  = (float)p.X.x;
        pa->dbg_y[p.fortran_number-1]  = (float)p.X.y;
        pa->dbg_z[p.fortran_number-1]  = (float)p.X.z();
        pa->dbg_px[p.fortran_number-1] = (float)p.pu;
        pa->dbg_py[p.fortran_number-1] = (float)p.pv;
        pa->dbg_pz[p.fortran_number-1] = (float)p.pw;
   	}

}

//template <template <class Particle,int dims> class Cell,int dims>
//__global__ void GPU_StepAllCells(Cell<Particle,dims>  **cells,
//		                         double mass,
//		                         double q_mass,
//		                         double *p_control,
//		                         int jmp,
//		                         int nt
//		                         )
//{
//	Cell<Particle,dims>  *c,*c0 = cells[0];
//	__shared__ extern CellDouble fd[9];
//	CellDouble *c_jx,*c_jy,*c_jz,*c_ex,*c_ey,*c_ez,*c_hx,*c_hy,*c_hz;
//	CurrentTensor t1,t2;
//	int pqr2;
//	Particle p;
//
//	c = cells[ c0->getGlobalCellNumber(blockIdx.x,blockIdx.y,blockIdx.z)];
//
//	c_ex = &(fd[0]);
//	c_ey = &(fd[1]);
//	c_ez = &(fd[2]);
//
//	c_hx = &(fd[3]);
//	c_hy = &(fd[4]);
//	c_hz = &(fd[5]);
//
//	c_jx = &(fd[6]);
//	c_jy = &(fd[7]);
//	c_jz = &(fd[8]);
//
//	int index  = threadIdx.x;
//
//
//	while(index < CellExtent*CellExtent*CellExtent)
//	{
//
//		copyCellDouble(c_ex,c->Ex,index,blockIdx);
//		copyCellDouble(c_ey,c->Ey,index,blockIdx);
//		copyCellDouble(c_ez,c->Ez,index,blockIdx);
//
//		copyCellDouble(c_hx,c->Hx,index,blockIdx);
//		copyCellDouble(c_hy,c->Hy,index,blockIdx);
//		copyCellDouble(c_hz,c->Hz,index,blockIdx);
//
//		copyCellDouble(c_jx,c->Jx,index,blockIdx);
//		copyCellDouble(c_jy,c->Jy,index,blockIdx);
//		copyCellDouble(c_jz,c->Jz,index,blockIdx);
//		index += blockDim.x;
//	}
//	__syncthreads();
//
//	index  = threadIdx.x;
//
//    while(index < c->number_of_particles)
//    {
//
//        c->Move (index,&pqr2,&t1,&t2,mass,q_mass,p_control,jmp,c_ex,c_ey,c_ez,c_hx,c_hy,c_hz);
//
//        writeCurrentComponent(c_jx,&(t1.Jx),&(t2.Jx),pqr2);
//        writeCurrentComponent(c_jy,&(t1.Jy),&(t2.Jy),pqr2);
//        writeCurrentComponent(c_jz,&(t1.Jz),&(t2.Jz),pqr2);
//
//        index += blockDim.x;
//    }
//    __syncthreads();
//    index  = threadIdx.x;
//
//	while(index < CellExtent*CellExtent*CellExtent)
//	{
//      	copyCellDouble(c->Jx,c_jx,index,blockIdx);
//    	copyCellDouble(c->Jy,c_jy,index,blockIdx);
//    	copyCellDouble(c->Jz,c_jz,index,blockIdx);
//
//    	index += blockDim.x;
//    }
//    c->busyParticleArray = 0;
//}




__device__ SingleNodeFunctionType d_GPU_SetAllCurrentsToZero_SingleNode = GPU_SetAllCurrentsToZero_SingleNode;
SingleNodeFunctionType h_GPU_SetAllCurrentsToZero_SingleNode;

__device__ SingleNodeFunctionType d_GPU_getCellEnergy_SingleNode = GPU_getCellEnergy_SingleNode;
SingleNodeFunctionType h_GPU_getCellEnergy_SingleNode;

__device__ SingleNodeFunctionType d_GPU_SetFieldsToCells_SingleNode = GPU_SetFieldsToCells_SingleNode;
SingleNodeFunctionType h_GPU_SetFieldsToCells_SingleNode;

__device__ SingleNodeFunctionType d_GPU_WriteAllCurrents_SingleNode = GPU_WriteAllCurrents_SingleNode;
SingleNodeFunctionType h_GPU_WriteAllCurrents_SingleNode;

__device__ SingleNodeFunctionType d_GPU_GetCellNumbers_SingleNode = GPU_GetCellNumbers_SingleNode;
SingleNodeFunctionType h_GPU_GetCellNumbers_SingleNode;

__device__ SingleNodeFunctionType d_GPU_StepAllCells_SingleNode = GPU_StepAllCells_SingleNode;
SingleNodeFunctionType h_GPU_StepAllCells_SingleNode;

__device__ SingleNodeFunctionType d_GPU_Diagnose_SingleNode = GPU_Diagnose_SingleNode;
SingleNodeFunctionType h_GPU_Diagnose_SingleNode;

__device__ SingleNodeFunctionType d_GPU_eme_SingleNode = GPU_eme_SingleNode;
SingleNodeFunctionType h_GPU_eme_SingleNode;

__device__ SingleNodeFunctionType d_GPU_periodic_SingleNode = GPU_periodic_SingleNode;
SingleNodeFunctionType h_GPU_periodic_SingleNode;

__device__ SingleNodeFunctionType d_GPU_MakeDepartureLists_SingleNode = GPU_MakeDepartureLists_SingleNode;
SingleNodeFunctionType h_GPU_MakeDepartureLists_SingleNode;

__device__ SingleNodeFunctionType d_ArrangeFlights = ArrangeFlights;
SingleNodeFunctionType h_ArrangeFlights;

//GPU_GetFlownBeamNumber_SingleNode
__device__ SingleNodeFunctionType d_GPU_GetFlownBeamNumber_SingleNode = GPU_GetFlownBeamNumber_SingleNode;
SingleNodeFunctionType h_GPU_GetFlownBeamNumber_SingleNode;


__device__ SingleNodeFunctionType d_AddParticlesToCells = AddParticlesToCells;
SingleNodeFunctionType h_AddParticlesToCells;

__device__ SingleNodeFunctionType d_CurrentPeriodic = CurrentPeriodic;
SingleNodeFunctionType h_CurrentPeriodic;

__device__ SingleNodeFunctionType d_emh2 = emh2;
SingleNodeFunctionType h_emh2;

__device__ SingleNodeFunctionType d_add = GPU_add;
SingleNodeFunctionType h_add;

__device__ SingleNodeFunctionType d_emh1 = MagneticField_SingleNode;
SingleNodeFunctionType h_emh1;


#ifdef __CUDACC__
template <template <class Particle,int dims> class Cell,int dims>
__global__ void GPU_Universal_Kernel(Cell<Particle,dims>  **cells,
		KernelParams *params,SingleNodeFunctionType snf
		)
{
		unsigned int nx = blockIdx.x;
		unsigned int ny = blockIdx.y;
		unsigned int nz = blockIdx.z;
		unsigned int th_nx = threadIdx.x;
		unsigned int th_ny = threadIdx.y;
		unsigned int th_nz = threadIdx.z;

//		SingleNodeFunctionType sn = GPU_SetAllCurrentsToZero_SingleNode;

//		printf("launch %d %d %d %d %d %d\n",nx,ny,nz,th_nx,th_ny,th_nz);
		snf(cells,params,nx,ny,nz,th_nx,th_ny,th_nz);
//		GPU_SetAllCurrentsToZero_SingleNode(cells,nx,ny,nz);
}
#endif







#ifdef __CUDACC__
template <template <class Particle,int dims> class Cell,int dims>
__global__ void GPU_SetFieldsToCells(Cell<Particle,dims>  **cells,
        double *Ex,double *Ey,double *Ez,
        double *Hx,double *Hy,double *Hz
		)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	Cell<Particle,dims>  *c,*c0 = cells[0],nc;

	c = cells[ c0->getGlobalCellNumber(nx,ny,nz)];

	nc = *c;

	nc.readFieldsFromArrays(Ex,Ey,Ez,Hx,Hy,Hz,threadIdx);
}
#endif

__host__ __device__
double CheckArraySize(double* a, double* dbg_a,int size)
	{
	    int wrong = 0;
#ifdef CHECK_ARRAY_SIZE_DEBUG_PRINTS
	    printf("begin array checking1=============================\n");
#endif
	    for(int n = 0;n < size;n++)
	    {
	        if(fabs(a[n] - dbg_a[n]) > SIZE_TOLERANCE)
		{

#ifdef CHECK_ARRAY_SIZE_DEBUG_PRINTS
		   printf("n %5d %15.5e dbg %15.5e diff %15.5e wrong %10d \n",
				   n,a[n],dbg_a[n],fabs(a[n] - dbg_a[n]),wrong++);
#endif
		}
	    }
#ifdef CHECK_ARRAY_SIZE_DEBUG_PRINTS
	    printf("  end array checking=============================\n");
#endif

	    return (1.0-((double)wrong/(size)));
	}

#ifdef __CUDACC__
template <template <class Particle,int dims> class Cell,int dims>
__global__ void GPU_WriteAllCurrents(Cell<Particle,dims>  **cells,int n0,
		      double *jx,double *jy,double *jz,double *rho)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;

	Cell<Particle,dims>  *c,*c0 = cells[0],nc;
	__shared__ extern CellDouble fd[9];

	c = cells[ c0->getGlobalCellNumber(nx,ny,nz)];

	 nc = *c;

	             int i1,l1,k1;
	        	 i1 = threadIdx.x;
	        	 l1 = threadIdx.y;
	        	 k1 = threadIdx.z;
    	         int n = nc.getFortranCellNumber(nc.i+i1-1,nc.l+l1-1,nc.k+k1-1);

    	         if (n < 0 ) n = -n;
        		 double t,t_x,t_y;//,jx_p,jy_p;
		         //int i_f,l_f;//k_f;
		         t_x = nc.Jx->M[i1][l1][k1];
		         int3 i3 = nc.getCellTripletNumber(n);


		         MultiThreadAdd(&(jx[n]),t_x);
		         t_y= nc.Jy->M[i1][l1][k1];
		         MultiThreadAdd(&(jy[n]),t_y);
		         t = nc.Jz->M[i1][l1][k1];
		         MultiThreadAdd(&(jz[n]),t);

}
#endif

#ifdef __CUDACC__
template <template <class Particle,int dims> class Cell,int dims>
__global__ void GPU_WriteControlSystem(Cell<Particle,dims>  **cells)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	Cell<Particle,dims>  *c,*c0 = cells[0],nc;
	__shared__ extern CellDouble fd[9];

	c = cells[ c0->getGlobalCellNumber(nx,ny,nz)];

	 nc = *c;

	 nc.SetControlSystemToParticles();

}
#endif



#ifdef __CUDACC__
template <template <class Particle,int dims> class Cell,int dims>
__global__ void GPU_CollectStrayParticles(Cell<Particle,dims>  **cells,int nt)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;

	int busy;
	Particle p;
	int n;
	Cell<Particle,dims>  *c,*c0 = cells[0],nc,*new_c;

	c = cells[ n = c0->getGlobalCellNumber(nx,ny,nz)];

	for(int i = 0;i < c->number_of_particles; i++)
	{
		c->readParticleFromSurfaceDevice(i,&p);
#ifdef STRAY_DEBUG_PRINTS
		    printf("STRAY-BASIC step %d cell %3d %d %d sort %d particle %d FORTRAN %5d X: %15.5e < %15.5e < %15.5e Y: %15.5e < %15.5e < %15.5e Z: %15.5e < %15.5e < %15.5e \n",
		    		nt,c->i,c->l,c->k,(int)p.sort,i,p.fortran_number,
		    		c->x0,p.x,c->x0+c->hx,
		    		c->y0,p.y,c->y0+c->hy,
		    		c->z0,p.z,c->z0+c->hz
		    		);
#endif
		if(!c->isPointInCell(p.GetX()))// || (p.fortran_number == 753) )//|| (p.fortran_number == 10572))
		{
#ifdef STRAY_DEBUG_PRINTS

   			    printf("STRAY-OUT step %3d cell %3d %d %d sort %d particle %d FORTRAN %5d X: %15.5e < %25.17e < %15.5e \n",
   			    		nt,c->i,c->l,c->k,(int)p.sort,i,p.fortran_number,c->x0,p.x,c->x0+c->hx);
#endif
            int new_n = c->getPointCell(p.GetX());
            new_c = cells[new_n];

            while (atomicCAS(&(c->busyParticleArray),0,1)) {}
               c->removeParticleFromSurfaceDevice(i,&p,&(c->number_of_particles));
               //c->busyParticleArray = 0;
              atomicExch(&(c->busyParticleArray),0u);
               i--;


               while (atomicCAS(&(new_c->busyParticleArray),0,1)) {}

              new_c->Insert(p);
#ifdef STRAY_DEBUG_PRINTS

   			    printf("STRAY-INSERT step %d %3d %d %d sort %d particle %d FORTRAN %5d X: %15.5e < %25.17e < %15.5e \n",
   			    		nt,
   			    		new_c->i,new_c->l,new_c->k,(int)p.sort,i,p.fortran_number,new_c->x0,p.x,new_c->x0+new_c->hx);
#endif
              atomicExch(&(new_c->busyParticleArray),0u);

		}
	}
	c->printCellParticles("STRAY-FINAL",nt);

}
#endif


__device__ void copyCellDouble(CellDouble *dst,CellDouble *src,unsigned int n)
{
	if(n < CellExtent*CellExtent*CellExtent)
	{
		double *d_dst,*d_src;

		d_dst = (double *)(dst->getMp());
		d_src = (double *)(src->getMp());

		d_dst[n] = d_src[n];
//		printf("copyCellDouble %d Jx %25.15e %25.15e \n",n,d_dst[n], d_src[n]);

	}
}

__device__ void copyCellDouble(CellDouble *dst,CellDouble *src,unsigned int n,uint3 qq)
{
	if(n < CellExtent*CellExtent*CellExtent)
	{
		double *d_dst,*d_src;

		d_dst = (double *)(dst->getMp());
		d_src = (double *)(src->getMp());

		d_dst[n] = d_src[n];
	}
}

#ifdef __CUDACC__
template <template <class Particle,int dims> class Cell,int dims>
__global__ void GPU_GetCellNumbers(Cell<Particle,dims>  **cells,
		                         int *numbers)
{
		Cell<Particle,dims>  *c;//,nc;
		c = cells[blockIdx.x];

		numbers[blockIdx.x] = (*c).number_of_particles;
}
#endif


template <template <class Particle,int dims> class Cell,int dims>
__global__ void GPU_ControlAllCellsCurrents(Cell<Particle,dims>  **cells,int n,int i,CellDouble *jx,CellDouble *jy,CellDouble *jz)
{
	Cell<Particle,dims>  *c,*c0 = cells[0],nc;
	__shared__ extern CellDouble fd[9];
	c = cells[ n ];

	nc = *c;

#ifdef GPU_CONTROL_ALL_CELLS_CURRENTS_PRINT
#endif


}

__host__ __device__
void emh2_Element(
		Cell<Particle,DIMENSIONS> *c,
		int i,int l,int k,
		double *Q,double *H)
{
	int n  = c->getGlobalCellNumber(i,l,k);

	H[n] += Q[n];
}

__host__ __device__
void add_Element(
		Cell<Particle,DIMENSIONS> *c,
		int i,int l,int k,
		double B0,double *H)
{

	if(i > 102)
	{
		int qq = 0;
	}

	int n  = c->getGlobalCellNumber(i,l,k);
//    printf("addE n %10d i %5d l %5d k %5d %15.5e \n",n,i,l,k,H[n]);
	H[n] += B0;
//    printf("bddE n %10d i %5d l %5d k %5d %15.5e \n",n,i,l,k,H[n]);
}

#ifdef __CUDACC__
template <template <class Particle,int dims> class Cell,int dims>
__global__
void GPU_emh2(
		 Cell<Particle,dims>  **cells,
		 int i_s,int l_s,int k_s,
							double *Q,double *H
		)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	Cell<Particle,dims>  *c0 = cells[0];

	emh2_Element(c0,i_s+nx,l_s+ny,k_s+nz,Q,H);
}
#endif




__host__ __device__
void emh1_Element(
		Cell<Particle,DIMENSIONS> *c,
		int i,int l,int k,
		double *Q,double *H,double *E1, double *E2,
		double c1,double c2,
		int dx1,int dy1,int dz1,int dx2,int dy2,int dz2)
{

    int n  = c->getGlobalCellNumber(i,l,k);
	int n1 = c->getGlobalCellNumber(i+dx1,l+dy1,k+dz1);
	int n2 = c->getGlobalCellNumber(i+dx2,l+dy2,k+dz2);

	double e1_n1 = E1[n1];
	double e1_n  = E1[n];
	double e2_n2 = E2[n2];
	double e2_n  = E2[n];

	double t  = 0.5*(c1*(e1_n1 - e1_n)- c2*(e2_n2 - e2_n));
    Q[n] = t;
    H[n] += Q[n];
}

#ifdef __CUDACC__
template <template <class Particle,int dims> class Cell,int dims>
__global__
void GPU_emh1(
		 Cell<Particle,dims>  **cells,
				            int i_s,int l_s,int k_s,
							double *Q,double *H,double *E1, double *E2,
							double c1,double c2,
							int dx1,int dy1,int dz1,int dx2,int dy2,int dz2
		)
{
	unsigned int nx = blockIdx.x;
	unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	Cell<Particle,dims>  *c0 = cells[0];

	emh1_Element(c0,i_s+nx,l_s+ny,k_s+nz,Q,H,E1,E2,c1,c2,dx1,dy1,dz1,dx2,dy2,dz2);
}
#endif

__host__ __device__
	void emeElement(Cell<Particle,DIMENSIONS> *c,int i,int l,int k,double *E,double *H1, double *H2,
			double *J,double c1,double c2, double tau,
			int dx1,int dy1,int dz1,int dx2,int dy2,int dz2
			)
	{
	   int n  = c->getGlobalCellNumber(i,l,k);
	  int n1 = c->getGlobalCellNumber(i+dx1,l+dy1,k+dz1);
	  int n2 = c->getGlobalCellNumber(i+dx2,l+dy2,k+dz2);

	  E[n] += c1*(H1[n] - H1[n1]) - c2*(H2[n] - H2[n2]) - tau*J[n];
	}

__host__ __device__
void periodicElement(Cell<Particle,DIMENSIONS> *c,int i,int k,double *E,int dir, int to,int from)
{
    int n   = c->getGlobalBoundaryCellNumber(i,k,dir,to);
	int n1  = c->getGlobalBoundaryCellNumber(i,k,dir,from);

	if(dir == 0)
	{
		E[n1] = 0.0;
	}
    E[n]    = E[n1];

}

#ifdef __CUDACC__
template <template <class Particle,int dims> class Cell,int dims>
__global__ void GPU_periodic(Cell<Particle,dims>  **cells,
                             int i_s,int k_s,
                             double *E,int dir, int to,int from)
{
	unsigned int nx = blockIdx.x;
	//unsigned int ny = blockIdx.y;
	unsigned int nz = blockIdx.z;
	Cell<Particle,dims>  *c0 = cells[0];

	periodicElement(c0,nx+i_s,nz+k_s,E, dir,to,from);
}
#endif


__host__ __device__
void periodicCurrentElement(Cell<Particle,DIMENSIONS> *c,int i,int k,double *E,int dir, int dirE,int N)
{
	if(dir == 0) return;

    int n1    = c->getGlobalBoundaryCellNumber(i,k,dir,1);
    int n_Nm1 = c->getGlobalBoundaryCellNumber(i,k,dir,N-1);
    if(dir != dirE)
    {
       E[n1] += E[n_Nm1];
    }
    if(dir != 1 || dirE != 1)
    {
       E[n_Nm1] =  E[n1];
    }

    int n_Nm2 = c->getGlobalBoundaryCellNumber(i,k,dir,N-2);
    int n0    = c->getGlobalBoundaryCellNumber(i,k,dir,0);

#ifdef PERIODIC_CURRENT_PRINTS
    printf("%e %e \n",E[n0],E[n_Nm2]);
#endif
    E[n0] += E[n_Nm2];
    E[n_Nm2] = E[n0];
}

#ifdef __CUDACC__
template <template <class Particle,int dims> class Cell,int dims>
__global__ void GPU_CurrentPeriodic(Cell<Particle,dims>  **cells,double *E,int dirE, int dir,
                             int i_s,int k_s,int N)
{
	unsigned int nx = blockIdx.x;
	unsigned int nz = blockIdx.z;
	Cell<Particle,dims>  *c0 = cells[0];

	periodicCurrentElement(c0,nx+i_s,nz+k_s,E, dir,dirE,N);
}
#endif

#ifdef __CUDACC__
template <template <class Particle,int dims> class Cell,int dims>
__global__ void GPU_eme(

		            Cell<Particle,dims>  **cells,
		            int i_s,int l_s,int k_s,
					double *E,double *H1, double *H2,
					double *J,double c1,double c2, double tau,
					int dx1,int dy1,int dz1,int dx2,int dy2,int dz2
		)
{
	unsigned int nx = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int ny = blockIdx.y*blockDim.y + threadIdx.y;
	unsigned int nz = blockIdx.z*blockDim.z + threadIdx.z;
	Cell<Particle,dims>  *c0 = cells[0];




	emeElement(c0,i_s+nx,l_s+ny,k_s+nz,E,H1,H2,
			    	  		J,c1,c2,tau,
			    	  		dx1,dy1,dz1,dx2,dy2,dz2);
}
#endif

#ifdef __CUDACC__
template <template <class Particle,int dims> class Cell,int dims>
__global__ void copy_pointers(Cell<Particle,dims>  **cells,int *d_flags,double_pointer *d_pointers)
{
	Cell<Particle,dims>  *c = cells[blockIdx.x];

	c->flag_wrong_current_cell = d_flags[blockIdx.x];
	c->d_wrong_current_particle_attributes = d_pointers[blockIdx.x];

}
#endif



template <template <class Particle,int dims> class Cell,int dims >
class GPUPlasma
{
public:
   int  flag3D;
   int beam_plasma,nt_start_from_file,start_phase;
   int total_steps,minor_steps;

   double plasma_dim_y,plasma_dim_z;

   KernelParams params,*d_params;
   int N,lp;

   char  unique_variant_name[200];
   float *diag_x,*diag_y,*diag_pu,*diag_pv,*diag_pw,*diag_m,*diag_q_m,*diag_z;

   double beam_length,beam_width;
   double lx,ly,lz,Tb,rimp,rbd,tex0,tey0,tez0,tol;

   double flown_beam_particles;

   ParticleArraysGroup initial;
   ParticleFloatArraysGroup diagnostics,*d_diagnostics,host_copy_d_diagnostics;

   GPUCell<Particle,dims> **h_CellArray,**d_CellArray;
   GPUCell<Particle,dims> **cp;
//   thrust::device_vector<Cell<Particle> > *d_AllCells;
   double *d_Ex,*d_Ey,*d_Ez,*d_Hx,*d_Hy,*d_Hz,*d_Jx,*d_Jy,*d_Jz,*d_Rho,*d_npJx,*d_npJy,*d_npJz;
   double *d_Qx,*d_Qy,*d_Qz;
   double *dbg_x,*dbg_y,*dbg_z,*dbg_px,*dbg_py,*dbg_pz;
   int total_particles;

   int *h_controlParticleNumberArray;

   int  jx_wrong_points_number;
   int3 *jx_wrong_points,*d_jx_wrong_points;

   double *ctrlParticles,*d_ctrlParticles,*check_ctrlParticles;

   ParticleTarget<Particle,DIMENSIONS> *tgt_right;

   //external magnetic field
   double Bx0,d_Bx0;

   int jmp,size_ctrlParticles;
   double ami,amb,amf;
   int real_number_of_particle[3];
   FILE *f_prec_report;

   //diagnostics
   int component_total,sorts;
//   float *diagnose,*d_diagnose;



   int CPU_field;

void SetNumberOfSteps(int ts,int ms,int ph)
{
	total_steps = ts;
	minor_steps = ms;
	start_phase = ph;
}

void SetPlasmaSize(double pl_y,double pl_z)
{
	plasma_dim_y = pl_y;
	plasma_dim_z = pl_z;
}

void change_work_dir()
{
#ifdef SEPARATE_DIR

	struct stat st = {0};

	if (stat(unique_variant_name, &st) != -1)
	{
//		puts("before chdir");
	    chdir(unique_variant_name);
	    printf("%d changed dir to %s\n",getRank(),unique_variant_name);
	}

#endif

}

void get_variant_name(char *name)
{
	char hostname[1000],day[1000];
	string s;
	
	char variant_name[200];

#ifndef DEBUG_PATH

    if(getRank() == 0)
    {
    	gethostname(hostname,1000);
	    s = hostname;
    	time_t t = time(NULL);
     	struct tm tm = *localtime(&t);
	    int sec = 3600*tm.tm_hour+60*tm.tm_min+tm.tm_sec;

	    sprintf(day,"%4d-%02d-%02d-%05d-on-",
			tm.tm_year-100+2000,
			tm.tm_mon+1,tm.tm_mday,sec);

	    strcat(day,hostname);
	    strcpy(variant_name,day);

	    mkdir(variant_name, 0700);
    }


//    printf("%d send dir name %s \n",getRank(),name);
    SendComputationName(variant_name);
//    printf("%d received dir name %s \n",getRank(),name);
    strcpy(unique_variant_name,variant_name);

#else
		strcpy(name,"DEBUG");
#endif
   if(nt_start_from_file < 0)
   {
      change_work_dir();
   }
}

void Initialize()
{
	get_variant_name(unique_variant_name);

	f_prec_report = fopen("control_points.dat","wt");
	fclose(f_prec_report);

	InitializeCPU();

    double cf = checkFields(Ex,Ey,Ez,Hx,Hy,Hz,
    		               Jx,Jy,Jz,Qx,Qy,Qz,
                           0,mesh.size2(),001,0);

    InitGPUParticles();
    InitGPUFields();

}



#ifdef __CUDACC__
__host__ __device__
#endif
int zdim(){return ((flag3D == 1)*mesh.dimz2() + (flag3D == 0));}




int Kernel_Launcher(
		Cell<Particle,dims>  **cells,KernelParams *params,
		unsigned int grid_size_x,unsigned int grid_size_y,unsigned int grid_size_z,
        unsigned int block_size_x,unsigned int block_size_y,unsigned int block_size_z,
        int shmem_size,
        SingleNodeFunctionType h_snf,string name
		)
{
	struct timeval tv1,tv2;
#ifdef __CUDACC__
    dim3 blocks(grid_size_x,grid_size_y,grid_size_z),threads(block_size_x,block_size_y,block_size_z);

//    SingleNodeFunctionType h_sn;

   // cudaMemcpyFromSymbol( &h_sn, d_set_to0, sizeof( SingleNodeFunctionType ) );
//    printf("launch %d %d %d %d %d %d %s \n",i,j,k,i1,j1,k1,name);

    
    gettimeofday(&tv1,NULL);
    nvtxRangePushA(name.c_str());
    GPU_Universal_Kernel<<<blocks,threads,shmem_size>>>(cells,params,h_snf);
    DeviceSynchronize();
    nvtxRangePop();
    gettimeofday(&tv2,NULL);

#else
    char hostname[1000];
    gethostname(hostname,1000);

#ifdef OMP_OUTPUT
    printf("function %s executed on %s \n",name.c_str(),hostname);
#endif

    gettimeofday(&tv1,NULL);

#ifdef OMP_THREADS
    omp_set_num_threads(OMP_NUM_THREADS);

#pragma omp parallel for
#endif

    for(int i = 0;i < grid_size_x;i++)
    {
//#pragma omp parallel for
    	for(int j = 0;j < grid_size_y;j++)
    	{
//#pragma omp parallel for
    		for(int k = 0;k < grid_size_z;k++)
    		{
#ifdef DETAILED_CURRENTS_WRITE
    		    			if(!strcmp(name.c_str(),"WriteAllCurrents"))
    		    			{
    		    				char ftag[100];
    		    				sprintf(ftag,"bWriteAllCurrents_%d_%d_%d_",i,j,k);

    		    				write_field_component(params->nt,d_Jx,ftag,"Jx",mesh.size2());
    		    			}
#endif
//#pragma omp parallel for
    		    for(int i1 = 0;i1 < block_size_x;i1++)
    		    {
//#pragma omp parallel for
    		    	for(int j1 = 0;j1 < block_size_y;j1++)
    		    	{
//#pragma omp parallel for
    		    		for(int k1 = 0;k1 < block_size_z;k1++)
    		    		{
#ifdef OMP_THREADS
#ifdef OMP_DETAILED_OUTPUT
    		    			printf("launch %s %d %d %d %d %d %d %s  thread %d total %d \n",name.c_str(),
    		    					i,j,k,i1,j1,k1,name.c_str(),
    		    					omp_get_thread_num(),
    		    					omp_get_num_threads()
    		    					);
#endif
#endif
    		    			h_snf(cells,params,i,j,k,i1,j1,k1);
    		    		}
    		    	}
    		    }
#ifdef DETAILED_CURRENTS_WRITE
    		    			if(!strcmp(name.c_str(),"WriteAllCurrents"))
    		    			{
    		    				char ftag[100];
    		    				sprintf(ftag,"aWriteAllCurrents_%d_%d_%d_",i,j,k);

    		    				write_field_component(params->nt,d_Jx,ftag,"Jx",mesh.size2());
    		    			}
#endif
    		}
    	}
    }
    gettimeofday(&tv2,NULL);
#endif

#ifdef OMP_OUTPUT
    printf("function %s executed in %e sec \n",name.c_str(),
    		  tv2.tv_sec-tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)*1e-6
    		);
#endif




    int err = getLastError();
    if(err != 0)
    {
    	printf("error %s(%d) after %s \n",getErrorString((error_t)err), err, name.c_str());
    	exit(0);
    }

    return (int)err;

}


void SetAllCurrentsToZero(
		Cell<Particle,dims>  **cells,
		unsigned int grid_size_x,unsigned int grid_size_y,unsigned int grid_size_z,
        unsigned int block_size_x,unsigned int block_size_y,unsigned int block_size_z
		)
{
//    dim3 blocks(grid_size_x,grid_size_y,grid_size_z),threads(block_size_x,block_size_y,block_size_z);

//    SingleNodeFunctionType h_sn;



//    cudaMemcpyFromSymbol( &h_sn, d_GPU_SetAllCurrentsToZero_SingleNode, sizeof( SingleNodeFunctionType ) );



    Kernel_Launcher(cells,&params,grid_size_x,grid_size_y,grid_size_z,
                         block_size_x,block_size_y,block_size_z,16000,h_GPU_SetAllCurrentsToZero_SingleNode,
                         "SetAllCurrentsToZero");



}

int CopySingleNodeFunctionPointers()
{
	MemoryAllocate((void **)&d_params,sizeof(KernelParams));

	COPY_FUNCTION_POINTER( &h_GPU_SetAllCurrentsToZero_SingleNode,d_GPU_SetAllCurrentsToZero_SingleNode);
	COPY_FUNCTION_POINTER( &h_GPU_getCellEnergy_SingleNode,d_GPU_getCellEnergy_SingleNode);
	COPY_FUNCTION_POINTER( &h_GPU_SetFieldsToCells_SingleNode,d_GPU_SetFieldsToCells_SingleNode);
	COPY_FUNCTION_POINTER( &h_GPU_WriteAllCurrents_SingleNode,d_GPU_WriteAllCurrents_SingleNode );
	COPY_FUNCTION_POINTER( &h_GPU_GetCellNumbers_SingleNode,d_GPU_GetCellNumbers_SingleNode );
	COPY_FUNCTION_POINTER( &h_GPU_StepAllCells_SingleNode,d_GPU_StepAllCells_SingleNode );
	COPY_FUNCTION_POINTER( &h_GPU_Diagnose_SingleNode,d_GPU_Diagnose_SingleNode );
	COPY_FUNCTION_POINTER( &h_GPU_eme_SingleNode,d_GPU_eme_SingleNode );
	COPY_FUNCTION_POINTER( &h_GPU_periodic_SingleNode,d_GPU_periodic_SingleNode );
	COPY_FUNCTION_POINTER( &h_GPU_MakeDepartureLists_SingleNode,d_GPU_MakeDepartureLists_SingleNode );
	COPY_FUNCTION_POINTER( &h_ArrangeFlights,d_ArrangeFlights );
	COPY_FUNCTION_POINTER( &h_GPU_GetFlownBeamNumber_SingleNode,d_GPU_GetFlownBeamNumber_SingleNode );
	COPY_FUNCTION_POINTER( &h_AddParticlesToCells,d_AddParticlesToCells );
	COPY_FUNCTION_POINTER( &h_CurrentPeriodic,d_CurrentPeriodic );
	COPY_FUNCTION_POINTER( &h_emh2,d_emh2 );
	COPY_FUNCTION_POINTER( &h_emh1,d_emh1 );

	COPY_FUNCTION_POINTER( &h_add,d_add );

	return 0;
}


void InitGPUFields()
{
	CopySingleNodeFunctionPointers();

	MemoryAllocate((void**)&d_Ex,sizeof(double)*mesh.size2());
	MemoryAllocate((void**)&d_Ey,sizeof(double)*mesh.size2());
	MemoryAllocate((void**)&d_Ez,sizeof(double)*mesh.size2());

	MemoryAllocate((void**)&d_Hx,sizeof(double)*mesh.size2());
	MemoryAllocate((void**)&d_Hy,sizeof(double)*mesh.size2());
	MemoryAllocate((void**)&d_Hz,sizeof(double)*mesh.size2());

	MemoryAllocate((void**)&d_Jx,sizeof(double)*mesh.size2());
	MemoryAllocate((void**)&d_Jy,sizeof(double)*mesh.size2());
	MemoryAllocate((void**)&d_Jz,sizeof(double)*mesh.size2());

	MemoryAllocate((void**)&d_npJx,sizeof(double)*mesh.size2());
	MemoryAllocate((void**)&d_npJy,sizeof(double)*mesh.size2());
	MemoryAllocate((void**)&d_npJz,sizeof(double)*mesh.size2());

	MemoryAllocate((void**)&d_Qx,sizeof(double)*mesh.size2());
	MemoryAllocate((void**)&d_Qy,sizeof(double)*mesh.size2());
	MemoryAllocate((void**)&d_Qz,sizeof(double)*mesh.size2());

    copyFieldsToGPU();
}

void copyFieldsToGPU()
{
	int err;

    err = MemoryCopy(d_Ex,Ex,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
    {
    	printf("1copyFieldsToGPU err %d %s \n",err,getErrorString(err));
    	exit(0);
    }
    err = MemoryCopy(d_Ey,Ey,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
    {
     	printf("2copyFieldsToGPU err %d %s \n",err,getErrorString(err));
    	exit(0);
    }

    err = MemoryCopy(d_Ez,Ez,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("3copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_Hx,Hx,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("4copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }
    err = MemoryCopy(d_Hy,Hy,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("5copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }
    err = MemoryCopy(d_Hz,Hz,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("6copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_Jx,Jx,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("7copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }
    err = MemoryCopy(d_Jy,Jy,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("8copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_Jz,Jz,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("9copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_npJx,npJx,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("10copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_npJy,npJy,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("11copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_npJz,npJz,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("12copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_Qx,Qx,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("13copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_Qy,Qy,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("14copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }

    err = MemoryCopy(d_Qz,Qz,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("15copyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }
}

void InitGPUParticles()
{
	int size;
	GPUCell<Particle,dims> *d_c,*h_ctrl;
	GPUCell<Particle,dims> *n;
	GPUCell<Particle,dims> *h_copy,*h_c;
	double t;
//	dim3 dimGrid(mesh.x+2,mesh.y+2,mesh.dimz2()),dimBlockOne(1,1,1);



	 readControlFile(START_STEP_NUMBER);

	size = (*AllCells).size();
	
	memoryAllocationLog(sizeof(int)*size,__FILE__,__LINE__);

	h_controlParticleNumberArray = (int *)malloc(size*sizeof(int));

	 size_t m_free,m_total;
	 
	 memoryAllocationLog(sizeof(Cell<Particle,dims>)*2,__FILE__,__LINE__);

	h_ctrl = new Cell<Particle,dims>;
	n = new Cell<Particle,dims>;

	memoryAllocationLog(size*sizeof(Cell<Particle,dims> *)*2,__FILE__,__LINE__);
	
    h_CellArray = (Cell<Particle,dims> **)malloc(size*sizeof(Cell<Particle,dims> *));
    int err = (int)MemoryAllocate((void**)&d_CellArray,size*sizeof(Cell<Particle,dims> *));


    printf("%s : size = %d\n", __FILE__, size);
    for(int i = 0;i < size;i++)
    {
    	if(i == 4115)
    	{
    		int q = 0;
    	}
    	GPUCell<Particle,dims> c;
    	c = (*AllCells)[i];
//    	printf("A i %10d AllCells %p \n",i,AllCells);

    	h_controlParticleNumberArray[i] = c.number_of_particles;
    	/////////////////////////////////////////
    	*n = c;
//    	printf("B i %10d AllCells %p \n",i,AllCells);
#ifdef ATTRIBUTES_CHECK
    	c.SetControlSystem(jmp,d_ctrlParticles);
#endif

    	if(i == 4115)
    	    	{
    	    		int q = 0;

    	    	}
//    	printf("C i %10d AllCells %p \n",i,AllCells);

        d_c = c.copyCellToDevice(AllCells);

        double mfree,mtot;
        mfree = m_free;
        mtot  = m_total;
//    	printf("D i %10d AllCells %p \n",i,AllCells);

#ifdef COPY_CELL_PRINTS
        printf("cell %10d Device cell array allocated error %d %s memory: free %10.2f total %10.2f\n",i,err,GetErrorString(err),
        		                                                mfree/1024/1024/1024,mtot/1024/1024/1024);
        puts("");

	  dbgPrintGPUParticleAttribute(d_c,50,1," CO2DEV " );
	  puts("COPY----------------------------------");
#endif

#ifdef PARTICLE_PRINTS

        if(t < 1.0)
        {
        	t = c.compareToCell(*h_copy);
        }
#endif
        ////////////////////////////////////////.
//    	printf("E i %10d AllCells %p \n",i,AllCells);

        h_CellArray[i] = d_c;
        MemoryCopy(h_ctrl,d_c,sizeof(Cell<Particle,dims>),DEVICE_TO_HOST);
//    	printf("F i %10d AllCells %p \n",i,AllCells);

#ifdef InitGPUParticles_PRINTS
	    dbgPrintGPUParticleAttribute(d_c,50,1," CPY " );

       cudaPrintfInit();

        testKernel<<<1,1>>>(h_ctrl->d_ctrlParticles,h_ctrl->jmp);
        cudaPrintfDisplay(stdout, true);
        cudaPrintfEnd();

        printf("i %d l %d k n %d %d %e src %e num %d\n",h_ctrl->i,h_ctrl->l,h_ctrl->k,i,
        		c.ParticleArrayRead(0,7),c.number_of_particles
        		);
	printf("GPU cell %d ended ******************************************************\n",i);
#endif
    }


    err = MemoryCopy(d_CellArray,h_CellArray,size*sizeof(Cell<Particle,dims> *),HOST_TO_DEVICE);
    if(err != 0)
        {
         	printf("bGPU_WriteControlSystem err %d %s \n",err,getErrorString(err));
        	exit(0);
        }


#ifdef ATTRIBUTES_CHECK
    GPU_WriteControlSystem<<<dimGrid, dimBlockOne,16000>>>(d_CellArray);
    err = getLastError();
    cudaThreadSyncronize();
    TEST_ERROR(err);
#endif
	size = 0;

}


void copyCells(string where,int nt)
{
	static int first = 1;
	size_t m_free,m_total;
	int size = (*AllCells).size();

	struct sysinfo info;
	unsigned long c1,c2;

    if(first == 1)
    {
        memoryAllocationLog(size*sizeof(Cell<Particle,dims> *),__FILE__,__LINE__);
	
    	cp = (Cell<Particle,dims> **)malloc(size*sizeof(Cell<Particle,dims> *));
    }

	unsigned long m1,m2,delta,accum;
	memory_monitor("beforeCopyCells",nt);

	for(int i = 0;i < size;i++)
	{
		int err = GetDeviceMemory(&m_free,&m_total);
		sysinfo(&info);
		m1 = info.freeram;
	 	GPUCell<Particle,dims> c,*d_c,*z0;
	 	z0 = h_CellArray[i];
	 	if(first == 1)
	 	{
	       d_c = c.allocateCopyCellFromDevice();
     	   cp[i] = d_c;
	 	}
	 	else
	 	{
	 	   d_c = cp[i];
	 	}
	    c.copyCellFromDevice(z0,d_c,where,nt);
		m2 = info.freeram;

	    delta = m1-m2;
        accum += delta;
	}

	if(first == 1)
	{
		first = 0;
	}

	memory_monitor("afterCopyCells",nt);


}

void freeCellCopies(Cell<Particle,dims> **cp)
{
	int size = (*AllCells).size();

	for(int i = 0;i < size;i++)
	{
		GPUCell<Particle,dims> *d_c,c;

		d_c = cp[i];

		c.freeCopyCellFromDevice(d_c);

	}
	free(cp);
}

double compareCells(int nt)
{
	double t = 0.0,t1;
	struct sysinfo info;

	int size = (*AllCells).size();
	checkParticleNumbers(cp,-1);
	memory_monitor("compareCells",nt);

	for(int i = 0;i < size;i++)
	{
	 	Cell<Particle,dims> c = (*AllCells)[i];
	    t1 = c.compareToCell(*(cp[i]));

	    Particle p;

	    c.readParticleFromSurfaceDevice(0,&p);

	    if(t1 < 1.0)
	    {
	       	t1 = c.compareToCell(*(cp[i]));
	    }
	    if(isNan(t1))
	    {
	    	t1 = c.compareToCell(*(cp[i]));
	    }

	    t += t1;

	}
	memory_monitor("compareCells2",nt);

	return t/size;
}

double checkGPUArray(double *a,double *d_a,char *name,char *where,int nt)
{
	 static double *t;
	 static int f1 = 1;
	 char fname[1000];
	 double res;
	 FILE *f;
#ifndef CHECK_ARRAY_OUTPUT
   return 0.0;
#endif

	 sprintf(fname,"diff_%s_at_%s_nt%08.dat",name,where,nt);


	 if(f1 == 1)
	 {
	         memoryAllocationLog(sizeof(double)*mesh.size2(),__FILE__,__LINE__);
		 t = (double *)malloc(sizeof(double)*mesh.size2());
		 f1 = 0;
	 }
	 int err;
	 err = MemoryCopy(t,d_a,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
	 if(err != 0)
	         {
	          	printf("bCheckArray err %d %s \n",err,getErrorString(err));
	        	exit(0);
	         }

	 if((f = fopen(fname,"wt")) != NULL)
	 {
// res = CheckArray(a,t,f);
		 fclose(f);
	 }

	 return res;
}

double checkGPUArray(double *a,double *d_a)
{
	 static double *t;
	 static int f = 1;

	 if(f == 1)
	 {
	         memoryAllocationLog(sizeof(double)*mesh.size2(),__FILE__,__LINE__);
		 t = (double *)malloc(sizeof(double)*mesh.size2());
		 f = 0;
	 }
	 int err;
	 err = MemoryCopy(t,d_a,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
	 if(err != 0)
	         {
	          	printf("bCheckArray err %d %s \n",err,getErrorString(err));
	        	exit(0);
	         }

	 return CheckArray(a,t);

}


double compareFields(){return 1.0;}

double compareCPUtoGPU()
{
	return 1.0;
}

void StepAllCells()
{
//	dim3 dimGrid(mesh.x+2,mesh.y+2,mesh.dimz2()),
//			dimBlock(MAX_particles_per_cell/2,1,1);

    int err;
    err = getLastError();
    printf("Err: %d %s\n", err, getErrorString(err));
    DeviceSynchronize();
    err = getLastError();
    printf("Err: %d %s\n", err, getErrorString(err));
}


public:



void virtual emeGPUIterate(int i_s,int i_f,int l_s,int l_f,int k_s,int k_f,
			double *E,double *H1, double *H2,
			double *J,double c1,double c2, double tau,
			int dx1,int dy1,int dz1,int dx2,int dy2,int dz2)
{
	int err;
//	dim3 dimGrid(i_f-i_s+1,1,1),dimBlock(1,l_f-l_s+1,k_f-k_s+1);

	params.i_s = i_s;
	params.k_s = k_s;
	params.l_s = l_s;

	params.E = E;
	params.H1 = H1;
	params.H2 = H2;
	params.J = J;
	params.c1 = c1;
	params.c2 = c2;
	params.tau = tau;

	params.dx1 = dx1;
	params.dy1 = dy1;
	params.dz1 = dz1;

	params.dx2 = dx2;
	params.dy2 = dy2;
	params.dz2 = dz2;
	params.blockDim_x = 1;
	params.blockDim_y = l_f-l_s+1;
	params.blockDim_z = k_f-k_s+1;


	MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

    Kernel_Launcher(d_CellArray,d_params,
    		             i_f-i_s+1,1,1,
    		             1,l_f-l_s+1,k_f-k_s+1,
    		             16000,
                         h_GPU_eme_SingleNode,
                         "eme");

//    GPU_eme<<<dimGrid,dimBlock>>>(d_CellArray,i_s,l_s,k_s,
//    		                            E,H1,H2,
//    					    	  		J,c1,c2,tau,
//    					    	  		dx1,dy1,dz1,dx2,dy2,dz2
//    		);
    err = getLastError();
    ThreadSynchronize();
//    TEST_ERROR(err);
}

int virtual ElectricFieldTrace(string lname,int nt,
  double *E,double *H1,double *H2,double *J,double *dbg_E0,double *dbg_E,double *dbg_H1,double *dbg_H2,double *dbg_J,int dir,double c1,double c2,double tau)
  {
      int i_start,l_start,k_start,dx1,dy1,dz1,dx2,dy2,dz2;
      double *dbg_E_aper;
      char logname[100];

      sprintf(logname,"%s%03d.dat",lname.c_str(),nt);


#ifdef DEBUG_PLASMA_EFIELDS
      dbg_E_aper = (double *)malloc(sizeof(double)*mesh.size2());

      read3DarrayLog(logname, dbg_E_aper,50,8);

      memoryAllocationLog(size*sizeof(sizeof(double)*mesh.size2()*5,__FILE__,__LINE__);


      Jloc = (double *)malloc(sizeof(double)*mesh.size2());
      ldH1 = (double *)malloc(sizeof(double)*mesh.size2());
      ldH2 = (double *)malloc(sizeof(double)*mesh.size2());
      ldE  = (double *)malloc(sizeof(double)*mesh.size2());


      char dbg_fnex[100];

     read3DarrayLog(logname, ldE,50,0);
     read3DarrayLog(logname, ldH1,50,1);
     read3DarrayLog(logname, ldH2,50,3);
     read3DarrayLog(logname, Jloc,50,5);

if(CPU_field)
{
     CheckArray(E,ldE);
     puts("Ex prev");
     //exit(0);
     CheckArray(Jloc,dbg_J);
     CheckArray(J,dbg_J);
     CheckArray(J,Jloc);
     puts("Jx");
     //exit(0);
     CheckArray(H1,ldH1);
     puts("H1");
     //exit(0);
     CheckArray(H2,ldH2);
     puts("H2");
     //exit(0);
}
else
{
    checkGPUArray(ldE,E);
    puts("Ex prev");
    //exit(0);
    checkGPUArray(Jloc,J);
    puts("Jx");
    //exit(0);
    checkGPUArray(ldH1,H1);
    puts("H1");
    //exit(0);
    checkGPUArray(ldH2,H2);
    puts("H2");
}
#endif

      i_start = (dir == 0)*0 + (dir == 1)*1 + (dir == 2)*1;
      l_start = (dir == 0)*1 + (dir == 1)*0 + (dir == 2)*1;
      k_start = (dir == 0)*1 + (dir == 1)*1 + (dir == 2)*0;

      dx1 = (dir == 0)*0    + (dir == 1)*0    + (dir == 2)*(-1);
      dy1 = (dir == 0)*(-1) + (dir == 1)*0    + (dir == 2)*0;
      dz1 = (dir == 0)*0    + (dir == 1)*(-1) + (dir == 2)*0;

      dx2 = (dir == 0)*0    + (dir == 1)*(-1) + (dir == 2)*0;
      dy2 = (dir == 0)*0    + (dir == 1)*0    + (dir == 2)*(-1);
      dz2 = (dir == 0)*(-1) + (dir == 1)*0    + (dir == 2)*0;

      if(CPU_field)
      {
         emeIterate(i_start,mesh.x,l_start,mesh.y,k_start,mesh.z(),
    		                E,H1,H2,
      		    	  		J,c1,c2,tau,
      		    	  		dx1,dy1,dz1,dx2,dy2,dz2);
      }
      else
      {
         emeGPUIterate(i_start,mesh.x,l_start,mesh.y,k_start,mesh.z(),
    	      		                E,H1,H2,
    	        		    	  		J,c1,c2,tau,
    	        		    	  		dx1,dy1,dz1,dx2,dy2,dz2);
      }


    return 0;
  }

void  ComputeField_FirstHalfStep(
		   double *locEx,double *locEy,double *locEz,
		   int nt,
		   double *locHx,double *locHy,double *locHz,
		   double *loc_npJx,double *loc_npJy,double *loc_npJz,
		   double *locQx,double *locQy,double *locQz)
{
	 double t_check[15];

	 CPU_field = 0;
	 if(CPU_field == 0)
	 {
		 t_check[0] = checkControlMatrix("emh1",nt,"qx",d_Qx);
		 t_check[1] = checkControlMatrix("emh1",nt,"qy",d_Qy);
		 t_check[2] = checkControlMatrix("emh1",nt,"qz",d_Qz);

		 t_check[3] = checkControlMatrix("emh1",nt,"ex",d_Ex);
		 t_check[4] = checkControlMatrix("emh1",nt,"ey",d_Ey);
		 t_check[5] = checkControlMatrix("emh1",nt,"ez",d_Ez);

		 t_check[6] = checkControlMatrix("emh1",nt,"hx",d_Hx);
		 t_check[7] = checkControlMatrix("emh1",nt,"hy",d_Hy);
		 t_check[8] = checkControlMatrix("emh1",nt,"hz",d_Hz);

		 emh1(d_Qx,d_Qy,d_Qz,d_Hx,d_Hy,d_Hz,nt,d_Ex,d_Ey,d_Ez);

		 t_check[9] = checkControlMatrix("emj1",nt,"qx",d_Qx);
		 t_check[10] = checkControlMatrix("emj1",nt,"qy",d_Qy);
		 t_check[11] = checkControlMatrix("emj1",nt,"qz",d_Qz);

		 t_check[12] = checkControlMatrix("emj1",nt,"hx",d_Hx);
		 t_check[13] = checkControlMatrix("emj1",nt,"hy",d_Hy);
		 t_check[14] = checkControlMatrix("emj1",nt,"hz",d_Hz);
	 }


	 CPU_field = 1;

	 emh1(locQx,locQy,locQz,locHx,locHy,locHz,nt,locEx,locEy,locEz);

#ifdef CONTROL_POINT_CHECK
	 checkControlPoint(50,nt,0);
#endif
}

virtual void ComputeField_SecondHalfStep(
		   double *locEx,double *locEy,double *locEz,
		   int nt,
		   double *locHx,double *locHy,double *locHz,
		   double *loc_npJx,double *loc_npJy,double *loc_npJz,
		   double *locQx,double *locQy,double *locQz)
{
#ifdef CPU_DEBUG_RUN

#ifdef CONTROL_POINT_CHECK
     checkControlPoint(275,nt,0);
#endif
     SetPeriodicCurrents(nt);

#ifdef CONTROL_POINT_CHECK
	 checkControlPoint(400,nt,0);
#endif
	 emh2(locHx,locHy,locHz,nt,locQx,locQy,locQz);
#endif

	 CPU_field = 0;
	 if(CPU_field == 0)
	 {
		emh2(d_Hx,d_Hy,d_Hz,nt,d_Qx,d_Qy,d_Qz);
#ifdef CPU_DEBUG_RUN
	    checkFirstHalfstep_emh2_GPUMagneticFields(nt);
#endif
	 }
#ifdef CPU_DEBUG_RUN

#ifdef CONTROL_POINT_CHECK
	 checkControlPoint(500,nt,0);
#endif

#endif
	     CPU_field = 0;
			 eme(d_Ex,d_Ey,d_Ez,nt,d_Hx,d_Hy,d_Hz,d_Jx,d_Jy,d_Jz);
#ifdef CPU_DEBUG_RUN
			  checkGPUSecondHalfstepFields(nt);

		 CPU_field = 1;
#ifdef FINAL_CONTROL_POINT_CHECK
		 if(nt == TOTAL_STEPS)
		 {
		    checkControlPoint(600,nt,0);
		 }
#endif
#endif


}

void eme(double *locEx,double *locEy,double *locEz,
		   int nt,
		   double *locHx,double *locHy,double *locHz,
		   double *loc_npJx,double *loc_npJy,double *loc_npJz)
{
      Cell<Particle,dims> c = (*AllCells)[0];
      double hx = c.get_hx(),hy = c.get_hy(),hz = c.get_hz();
      double c11 = tau/hx,c21 = tau/hy,c31 = tau/hz;

      ElectricFieldTrace("exlg",nt,locEx,locHz,locHy,loc_npJx,
    		  dbgEx0,npEx,dbgHz,
    		  dbgHy,dbgJx,0,c21,c31,tau);
       PeriodicBoundaries(locEx,1,0,mesh.x,1,mesh.z(),mesh.y);
       PeriodicBoundaries(locEx,2,0,mesh.x,0,mesh.y+1,mesh.z());
      ElectricFieldTrace("eylg",nt,locEy,locHx,locHz,loc_npJy,
    		  dbgEy0,npEy,dbgHx,
    		  dbgHz,dbgJy,1,c31,c11,tau);

       PeriodicBoundaries(locEy,0,0,mesh.y,1,mesh.z(),mesh.x);
       PeriodicBoundaries(locEy,2,0,mesh.x+1,0,mesh.y,mesh.z());
       SinglePeriodicBoundary(locEy,1,0,mesh.x+1,0,mesh.z()+1,mesh.y);

      ElectricFieldTrace("ezlg",nt,locEz,locHy,locHx,loc_npJz,
    		  dbgEz0,npEz,dbgHy,
    		  dbgHx,dbgJz,2,c11,c21,tau);
      PeriodicBoundaries(locEz,0,1,mesh.y,0,mesh.z(),mesh.x);
      PeriodicBoundaries(locEz,1,0,mesh.x+1,0,mesh.z(),mesh.y);
}


virtual void emh1(
                  double *locQx,double *locQy,double *locQz,
                  double *locHx,double *locHy,double *locHz,
		            int nt,
		            double *locEx,double *locEy,double *locEz
		           )
{
    Cell<Particle,dims> c = (*AllCells)[0];
    double hx = c.get_hx(),hy = c.get_hy(),hz = c.get_hz();
    double c11 = tau/(hx),c21 = tau/(hy),c31 = tau/hz;

#ifdef DEBUG_PLASMA_STEP_FIELDS
    memoryAllocationLog(size*sizeof(sizeof(double)*mesh.size2()*6,__FILE__,__LINE__);
    
    d_locEx = (double *)malloc(sizeof(double)*mesh.size2());
    d_locEy = (double *)malloc(sizeof(double)*mesh.size2());
    d_locEz = (double *)malloc(sizeof(double)*mesh.size2());

    d_locHx = (double *)malloc(sizeof(double)*mesh.size2());
    d_locHy = (double *)malloc(sizeof(double)*mesh.size2());
    d_locHz = (double *)malloc(sizeof(double)*mesh.size2());

    readDebugArray("dnex",d_locEx,2*nt-2);
    readDebugArray("dney",d_locEy,2*nt-2);
    readDebugArray("dnez",d_locEz,2*nt-2);
    CheckArray(d_locEx,locEx);
    CheckArray(d_locEy,locEy);
    CheckArray(d_locEz,locEz);

    readDebugArray("dnqx",dbg_Qx,2*nt-1);
    readDebugArray("dnqy",dbg_Qy,2*nt-1);
    readDebugArray("dnqz",dbg_Qz,2*nt-1);

    readDebugArray("dnhx",dbgHx,2*nt);
    readDebugArray("dnhy",dbgHy,2*nt);
    readDebugArray("dnhz",dbgHz,2*nt);

    readDebugArray("dnhx",locHx,2*nt-1);
    readDebugArray("dnhy",locHy,2*nt-1);
    readDebugArray("dnhz",locHz,2*nt-1);

    readDebugArray("dnhx",d_locHx,2*nt-1);
    readDebugArray("dnhy",d_locHy,2*nt-1);
    readDebugArray("dnhz",d_locHz,2*nt-1);
    CheckArray(d_locHx,locHx);
    CheckArray(d_locHy,locHy);
    CheckArray(d_locHz,locHz);

#endif

     MagneticFieldTrace(c,"hxlg",nt,locQx,locHx,locEy,locEz,mesh.x+1,mesh.y,mesh.z(),c31,c21,0);
#ifdef DEBUG_PLASMA_STEP_FIELDS
    CheckArray(Qx,dbg_Qx);
    CheckArray(Hx,dbgHx);
#endif

    MagneticFieldTrace(c,"hylg",nt,locQy,locHy,locEz,locEx,mesh.x,mesh.y+1,mesh.z(),c11,c31,1);
#ifdef DEBUG_PLASMA_STEP_FIELDS
    CheckArray(Qy,dbg_Qy);
    CheckArray(Hy,dbgHy);
#endif

    MagneticFieldTrace(c,"hzlg",nt,locQz,locHz,locEx,locEy,mesh.x,mesh.y,mesh.z()+1,c21,c11,2);
#ifdef DEBUG_PLASMA_STEP_FIELDS
    CheckArray(Qz,dbg_Qz);
    CheckArray(Hz,dbgHz);
#endif
}
virtual void emh2(double *locHx,double *locHy,double *locHz,
		            int nt,
		            double *locQx,double *locQy,double *locQz)
{
    Cell<Particle,dims> c = (*AllCells)[0];

#ifdef DEBUG_PLASMA_STEP_FIELDS_EMH2

    readDebugArray("dnqx",dbg_Qx,nt,0);
    readDebugArray("dnqy",dbg_Qy,nt,0);
    readDebugArray("dnqz",dbg_Qz,nt,0);
    CheckArray(dbg_Qx,locQx);
    CheckArray(dbg_Qy,locQy);
    CheckArray(dbg_Qz,locQz);

    readDebugArray("dnhx",dbgHx,nt+1,0);
    readDebugArray("dnhy",dbgHy,nt+1,0);
    readDebugArray("dnhz",dbgHz,nt+1,0);


    CheckArray(dbgHx,locHx);
    CheckArray(dbgHy,locHy);
    CheckArray(dbgHz,locHz);

#endif

    SimpleMagneticFieldTrace(c,locQx,locHx,mesh.x+1,mesh.y,mesh.z());
    SimpleMagneticFieldTrace(c,locQy,locHy,mesh.x,mesh.y+1,mesh.z());
    SimpleMagneticFieldTrace(c,locQz,locHz,mesh.x,mesh.y,mesh.z()+1);

//    printGPUArray	(locHx,501,nt,"Bx");
   // if(nt == 1)
   // {
    	AddConstantMagneticField(c,locHx,mesh.x+1,mesh.y,mesh.z(),Bx0,nt);
  //  }
//    printGPUArray	(locHx,510,nt,"Bx");

}



	void Step(int nt)
	 {
		 char fn[100];
		 double cf = checkFields(d_Ex,d_Ey,d_Ez,d_Hx,d_Hy,d_Hz,
				 d_Jx,d_Jy,d_Jz,d_Qx,d_Qy,d_Qz,
		                         1,mesh.size2(),010,0);

		// printGPUParticle(14536,2);

 		 if(nt == START_STEP_NUMBER)
 		 {
 			 readControlPoint(NULL,fn,0,nt,0,1,Ex,Ey,Ez,Hx,Hy,Hz,Jx,Jy,Jz,Qx,Qy,Qz,
 		                           dbg_x,dbg_y,dbg_z,dbg_px,dbg_py,dbg_pz
				 );
 			copyFieldsToGPU();
 		 }

#ifdef CONTROL_POINT_CHECK
		checkControlPoint(0,nt,1);
#endif
		 if(flagCPUandGPUrun)
		 {
        	 memory_monitor("beforeComputeField_FirstHalfStep",nt);

             CPU_field = 0;
			 ComputeField_FirstHalfStep(
					  Ex,Ey,Ez,
					  nt,
					  Hx,Hy,Hz,
					  npJx,npJy,npJz,
					  Qx,Qy,Qz);

				memset(Jx,0,sizeof(double)*mesh.size2());
			    memset(Jy,0,sizeof(double)*mesh.size2());
			    memset(Jz,0,sizeof(double)*mesh.size2());

		memory_monitor("afterComputeField_FirstHalfStep",nt);


			  AssignCellsToArraysGPU();
			  memory_monitor("before_CellOrder_StepAllCells",nt);

			  double mass = -1.0/1836.0,q_mass = -1.0;
			  double cf250 = checkFields(d_Ex,d_Ey,d_Ez,d_Hx,d_Hy,d_Hz,
			  			  			 d_Jx,d_Jy,d_Jz,d_Qx,d_Qy,d_Qz,
			  			  		     1,mesh.size2(),250,nt);
//			  printf("b CellOrder rank %5d \n",getRank());
		      CellOrder_StepAllCells(nt,	mass,q_mass,1);
			  double cf260 = checkFields(d_Ex,d_Ey,d_Ez,d_Hx,d_Hy,d_Hz,
			  			  			 d_Jx,d_Jy,d_Jz,d_Qx,d_Qy,d_Qz,
			  			  		     1,mesh.size2(),260,nt);

		      memory_monitor("after_CellOrder_StepAllCells",nt);

#ifdef ATTRIBUTES_CHECK
		      checkParticleAttributes(nt);
#endif

#ifdef CONTROL_POINT_CHECK
			  checkControlPoint(270,nt,1);
#endif

	   memory_monitor("before270",nt);

#ifdef CPU_DEBUG_RUN
		      ComputeField_SecondHalfStep(
					  Ex,Ey,Ez,
					  nt,
					  Hx,Hy,Hz,
					  npJx,npJy,npJz,
					  Qx,Qy,Qz);
#endif

	   memory_monitor("after_ComputeField_SecondHalfStep",nt);

		    CPU_field = 0;
		 }

//		 Diagnose(nt);
		 double cf1 = checkFields(d_Ex,d_Ey,d_Ez,d_Hx,d_Hy,d_Hz,
				 d_Jx,d_Jy,d_Jz,d_Qx,d_Qy,d_Qz,
				                  1,mesh.size2(),300,0);

	 }
	virtual double getElectricEnergy()
	{
		int err;
//		dim3 dimGrid(mesh.x+2,mesh.y+2,mesh.dimz2()),dimGridOne(1,1,1),dimBlock(MAX_particles_per_cell/2,1,1),
//    		 dimBlockOne(1,1,1),dimBlockGrow(1,1,1),dimBlockExt(CellExtent,CellExtent,CellExtent);
		static int first = 1;
		static double *d_ee;
		double ee;


//		cudaMemset(d_ee,0,sizeof(double));
		params.d_ee = 0;
		params.d_Ex = d_Ex;
		params.d_Ey = d_Ey;
		params.d_Ez = d_Ez;
//		static int first  = 1;


//			first = 0;

		MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

	    Kernel_Launcher(d_CellArray,d_params,mesh.x+2,mesh.y+2,mesh.dimz2(),1,1,1,
	                         16000,h_GPU_getCellEnergy_SingleNode,
	                         "getCellEnergy");

//		GPU_getCellEnergy<<<dimGrid, dimBlockOne,16000>>>(d_CellArray,d_ee,d_Ex,d_Ey,d_Ez);

	    err = getLastError();
	    ThreadSynchronize();
//	    TEST_ERROR(err);

        MemoryCopy(&params,d_params,sizeof(KernelParams),DEVICE_TO_HOST);


        return params.d_ee;

	}

	int write_fields(int nt)
	{
		static int first = 1;
		static double *h_Ex,*h_Ey,*h_Ez,*h_Bx,*h_By,*h_Bz,*h_Jx,*h_Jy,*h_Jz;

		if(first == 1)
		{
			first = 0;

			MemoryAllocate((void **)&h_Ex,sizeof(double)*mesh.size2());
			MemoryAllocate((void **)&h_Ey,sizeof(double)*mesh.size2());
			MemoryAllocate((void **)&h_Ez,sizeof(double)*mesh.size2());

			MemoryAllocate((void **)&h_Bx,sizeof(double)*mesh.size2());
			MemoryAllocate((void **)&h_By,sizeof(double)*mesh.size2());
			MemoryAllocate((void **)&h_Bz,sizeof(double)*mesh.size2());

			MemoryAllocate((void **)&h_Jx,sizeof(double)*mesh.size2());
			MemoryAllocate((void **)&h_Jy,sizeof(double)*mesh.size2());
			MemoryAllocate((void **)&h_Jz,sizeof(double)*mesh.size2());
		}

		MemoryCopy(h_Ex,d_Ex,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
		MemoryCopy(h_Ey,d_Ey,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
		MemoryCopy(h_Ez,d_Ez,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);

		MemoryCopy(h_Bx,d_Hx,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
		MemoryCopy(h_By,d_Hy,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
		MemoryCopy(h_Bz,d_Hz,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);

		MemoryCopy(h_Jx,d_Jx,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
		MemoryCopy(h_Jy,d_Jy,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
		MemoryCopy(h_Jz,d_Jz,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);

		 Cell<Particle,dims> c = (*AllCells)[0];
		 double hx = c.get_hx(),hy = c.get_hy(),hz = c.get_hz();
		 char fname[100];
		 FILE *f;



		 sprintf(fname,"fields_%s_nt%08d.dat",unique_variant_name,nt);

	 if(getRank() == 0)
	 {
		 if((f = fopen(fname,"wt")) == NULL) return 1;

		 for(int i = 0;i < c.mesh.x;i++)
		 {
			 for(int l = 0;l < c.mesh.y;l++)
			 {
				 for(int k = 0;k < c.mesh.z();k++)
				 {
					 int n  = c.getGlobalCellNumber(i,l,k);

					 fprintf(f,"%10.3e %10.3e %10.3e %25.15e %25.15e %25.15e %25.15e %25.15e %25.15e %25.15e %25.15e %25.15e \n",
							 i*hx,l*hy,k*hz,
							 h_Ex[n],h_Ey[n],h_Ez[n],
							 h_Bx[n],h_By[n],h_Bz[n],
							 h_Jx[n],h_Jy[n],h_Jz[n]
							 );
				 }
			 }
		 }
		 fclose(f);
	 }

		 return 0;
	}

	void Diagnose(int nt)
	{
		double ee;
		static FILE *f;
		static int first = 1;

		if(first == 1)
		{
			f = fopen("energy.dat","wt");
			first = 0;
		}
		else
		{
			f = fopen("energy.dat","at");

		}

		if(nt % minor_steps == 0)
		{

		 params.component_total = component_total;
		 params.sorts = sorts;
		 params.d_part_diag = d_diagnostics;

		 MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

//		 GPU_Diagnose_SingleNode(d_CellArray,d_params,0,0,0, 0,0,0);
//		 GPU_Diagnose_SingleNode(d_CellArray,d_params,70,0,0, 0,0,0);


//		 Kernel_Launcher(d_CellArray,d_params,mesh.x+2,mesh.y+2,mesh.dimz2(),
//				    		1,1,1,
//				    		16000,h_GPU_Diagnose_SingleNode,
//		//		    		PARTICLE_BLOCK_X_SIZE,1,1,16000,h_GPU_StepAllCells_SingleNode,
//				            "Diagnose");
//
////		 MemoryCopy(diagnose,d_diagnose,sizeof(float)*sorts*component_total*8,DEVICE_TO_HOST);
//		 CopyParticleDiagnosticPointersFromDevice(diagnostics,
//				 host_copy_d_diagnostics);

//		 TransformParticlesToLists(SORTS,
//		        initial,
//		        diagnose,
//		        diagnostics);
//		 write_particles(&diagnostics,&initial,unique_variant_name,nt);
		 write_fields(nt);





//		 Particle2D_Distributions(sorts,unique_variant_name,nt,diagnostics,xmax.x,xmax.y);
		}

//		 tgt_right->Diagnose(nt);




        ee = getElectricEnergy();

        	fprintf(f,"%10d %25.15e \n",nt,ee);

        fclose(f);

	}
	virtual ~GPUPlasma(){
		}

      Point<double,DIMENSIONS,0> xmax,beam_max;
	  Point<int,DIMENSIONS,1> mesh;

//	  int mesh.x,mesh.y,mesh.z();

	  int n_per_cell;

	  int meh;

	  int magf;

	  double ion_q_m,tau;

//	  double Lx,Ly,Lz;

	  double ni;

	  double *Qx,*Qy,*Qz,*dbg_Qx,*dbg_Qy,*dbg_Qz;

	  double *Ex,*Ey,*Ez,*Hx,*Hy,*Hz,*Jx,*Jy,*Jz,*Rho,*npJx,*npJy,*npJz;
	  double *dbgEx,*dbgEy,*dbgEz,*dbgHx,*dbgHy,*dbgHz,*dbgJx,*dbgJy,*dbgJz;
	  double *dbgEx0,*dbgEy0,*dbgEz0;
	  double *npEx,*npEy,*npEz;



	  vector<Cell<Particle,dims> > *AllCells;

	  int getBoundaryLimit(int dir){return ((dir == 0)*mesh.x  + (dir == 1)*mesh.y + (dir == 2)*mesh.z() + 2);}

	  virtual void Alloc()
	  {

		  AllCells = new vector<Cell<Particle,dims> >;

		  int s2 = mesh.size2();
		  
		  memoryAllocationLog(mesh.size2()*sizeof(double)*34,__FILE__,__LINE__);

	     Ex  = new double[mesh.size2()];
	     Ey  = new double[mesh.size2()];
	     Ez  = new double[mesh.size2()];
	     Hx  = new double[mesh.size2()];
	     Hy  = new double[mesh.size2()];
	     Hz  = new double[mesh.size2()];
	     Jx  = new double[mesh.size2()];
	     Jy  = new double[mesh.size2()];
	     Jz  = new double[mesh.size2()];
	     Rho = new double[mesh.size2()];

	     npJx  = new double[mesh.size2()];
	     npJy  = new double[mesh.size2()];
	     npJz  = new double[mesh.size2()];

	     npEx  = new double[mesh.size2()];
	     npEy  = new double[mesh.size2()];
	     npEz  = new double[mesh.size2()];

	     Qx  = new double[mesh.size2()];
	     Qy  = new double[mesh.size2()];
	     Qz  = new double[mesh.size2()];

	#ifdef DEBUG_PLASMA

	     dbgEx  = new double[mesh.size2()];
	     dbgEy  = new double[mesh.size2()];
	     dbgEz  = new double[mesh.size2()];
	     dbgEx0  = new double[mesh.size2()];
	     dbgEy0  = new double[mesh.size2()];
	     dbgEz0  = new double[mesh.size2()];

	     dbgHx  = new double[mesh.size2()];
	     dbgHy  = new double[mesh.size2()];
	     dbgHz  = new double[mesh.size2()];
	     dbgJx  = new double[mesh.size2()];
	     dbgJy  = new double[mesh.size2()];
	     dbgJz  = new double[mesh.size2()];

	     dbg_Qx  = new double[mesh.size2()];
	     dbg_Qy  = new double[mesh.size2()];
	     dbg_Qz  = new double[mesh.size2()];
	#endif
	  }

	  virtual void InitFields()
	  {
	     for(int i = 0;i < mesh.size2();i++)
	     {
	         Ex[i] = 0.0;
	         Ey[i] = 0.0;
	         Ez[i] = 0.0;
	         Hx[i] = 0.0;
	         Hy[i] = 0.0;
	         Hz[i] = 0.0;

	         dbgEx[i] = 0.0;
	         dbgEy[i] = 0.0;
	         dbgEz[i] = 0.0;
	         dbgHx[i] = 0.0;
	         dbgHy[i] = 0.0;
	         dbgHz[i] = 0.0;
	     }
	  }

	  virtual void InitCells()
	  {
	     for(int i = 0;i < mesh.x+2;i++)
	     {
	         for(int l = 0;l < mesh.y+2;l++)
		 {
		     for(int k = 0;k < mesh.dimz2();k++)
		     {
	                 Cell<Particle,dims> * c = new Cell<Particle,dims>(
	                		 i,l,k,
	                		 xmax.x,xmax.y,xmax.z(),
	                		 mesh.x,mesh.y,mesh.z(),tau);
	                 c->Init();
	                 c->beam_boundary = beam_length;
	                 c->tgt           = tgt_right;

			         (*AllCells).push_back(*c);
#ifdef INIT_CELLS_DEBUG_PRINT
	                printf("%5d %5d %5d size %d \n",i,l,k,(*AllCells).size());
#endif
		     }

		 }

	     }
	  }

	  virtual void InitCurrents()
	  {
	     for(int i = 0;i < mesh.size2();i++)
	     {
	         Jx[i]  = 0.0;
	         Jy[i]  = 0.0;
	         Jz[i]  = 0.0;
	         Rho[i] = 0.0;

	         dbgJx[i]  = 0.0;
	         dbgJy[i]  = 0.0;
	         dbgJz[i]  = 0.0;

	     }
	  }

	  void InitCurrents(char *fnjx,char *fnjy,char *fnjz,
	                    char *dbg_fnjx,char *dbg_fnjy,char *dbg_fnjz,
	                    char *np_fnjx,char *np_fnjy,char *np_fnjz,
			            int dbg)
	  {
	     read3Darray(np_fnjx, npJx);
	     read3Darray(np_fnjy, npJy);
	     read3Darray(np_fnjz, npJz);

	     if(dbg == 0)
	     {
	        read3Darray(fnjx, Jx);
	        read3Darray(fnjy, Jy);
	        read3Darray(fnjz, Jz);
	     }
	#ifdef DEBUG_PLASMA
	     read3Darray(dbg_fnjx, dbgJx);
	     read3Darray(dbg_fnjy, dbgJy);
	     read3Darray(dbg_fnjz, dbgJz);

	#endif
	   }

	  void InitFields(char *fnex,char *fney,char *fnez,
			          char *fnhx,char *fnhy,char *fnhz,
			          char *dbg_fnex,char *dbg_fney,char *dbg_fnez,
			          char *dbg_0fnex,char *dbg_0fney,char *dbg_0fnez,
			          char *np_ex,char *np_ey,char *np_ez,
			          char *dbg_fnhx,char *dbg_fnhy,char *dbg_fnhz)
	  {
	     InitFields();

	     read3Darray(fnex, Ex);
	     read3Darray(fney, Ey);
	     read3Darray(fnez, Ez);
	     read3Darray(fnhx, Hx);
	     read3Darray(fnhy, Hy);
	     read3Darray(fnhz, Hz);

	#ifdef DEBUG_PLASMA
	     read3Darray(dbg_fnex, dbgEx);
	     read3Darray(dbg_fney, dbgEy);
	     read3Darray(dbg_fnez, dbgEz);

	     read3Darray(dbg_0fnex, dbgEx0);
	     read3Darray(dbg_0fney, dbgEy0);
	     read3Darray(dbg_0fnez, dbgEz0);

	     read3Darray(dbg_fnhx, dbgHx);
	     read3Darray(dbg_fnhy, dbgHy);
	     read3Darray(dbg_fnhz, dbgHz);

	     read3DarrayLog(np_ex, npEx,50,8);
	     read3DarrayLog(np_ey, npEy,50,8);
	     read3DarrayLog(np_ez, npEz,50,8);
	#endif
	  }


	  virtual void InitParticles(vector<Particle> & vp)
	  {
//	     InitIonParticles(n_per_cell,ion_q_m,vp);
	  }

	  virtual void InitParticles(char *fname,vector<Particle>& vp)
	  {
	     FILE *f;
	     char str[1000];
	     double x,y,z,px,py,pz,q_m,m;
	     int n = 0;

	     if((f = fopen(fname,"rt")) == NULL) return;

	     while(fgets(str,1000,f) != NULL)
	     {
	          x   = atof(str);
	          y   = atof(str + 25);
	          z   = atof(str + 50);
	          px  = atof(str + 75);
	          py  = atof(str + 100);
	          pz  = atof(str + 125);
	          m   = fabs(atof(str + 150));
	          q_m = atof(str + 175);
	#undef GPU_PARTICLE
		  memoryAllocationLog(sizeof(Particle),__FILE__,__LINE__);
		  Particle *p = new Particle(x,y,z,px,py,pz,m,q_m);
		  p->fortran_number = ++n;
		  vp.push_back(*p);
	#define GPU_PARTICLE

	     }

         memoryAllocationLog(sizeof(double)*vp.size()*6,__FILE__,__LINE__);
         dbg_x = (double *)malloc(sizeof(double)*vp.size());
         dbg_y = (double *)malloc(sizeof(double)*vp.size());
         dbg_z = (double *)malloc(sizeof(double)*vp.size());
         dbg_px = (double *)malloc(sizeof(double)*vp.size());
         dbg_py = (double *)malloc(sizeof(double)*vp.size());
         dbg_pz = (double *)malloc(sizeof(double)*vp.size());

         total_particles = vp.size();

	     magf = 1;
	  }


	  void debugPrintParticleCharacteristicArray(double *p_ch,int np,int nt,char *name,int sort)
	  {
		   char fname[200];
		   FILE *f;

#ifndef PRINT_PARTICLE_INITIALS
		   return;

#else
		   sprintf(fname,"particle_init_%s_%05d_sort%02d.dat",name,nt,sort);

		   if((f = fopen(fname,"wt")) == NULL) return;

		   for (int i = 0;i < np;i++)
		   {
			   fprintf(f,"%10d %10d %25.16e \n",i,i+1,p_ch[i]);
		   }
		   fclose(f);
#endif
	  }

//      virtual int readBinaryParticleArraysOneSort(
//    		  FILE *f,
//    		  double **dbg_x,
//    		  double **dbg_y,
//    		  double **dbg_z,
//    		  double **dbg_px,
//    		  double **dbg_py,
//    		  double **dbg_pz,
//    		  double *qq_m,
//    		  double *mm,
//    		  int nt,
//    		  int sort
//    		  )
//      {
//		     double q_m,tp,m;
//		     int t;
//		     Cell<Particle,dims> c0 = (*AllCells)[0];
//		     int total_particles;
//		     int err;
//
//		     if((err = ferror(f)) != 0)
//		     {
//		     	 return err ;
//		     }
//
//		     fread(&t,sizeof(int),1,f);
//		     if((err = ferror(f)) != 0)
//		    	 {
//		    	 	 return err ;
//		    	 }
//		     fread(&tp,sizeof(double),1,f);
//		     if((err = ferror(f)) != 0)
//		    	 {
//		    	 	 return err ;
//		    	 }
//
//		     total_particles = (int)tp;
//		     fread(&q_m,sizeof(double),1,f);
//		     if((err = ferror(f)) != 0)
//		    	 {
//		    	 	 return err ;
//		    	 }
//
//		     fread(&m,sizeof(double),1,f);
//		     if((err = ferror(f)) != 0)
//		    	 {
//		    	 	 return err ;
//		    	 }
//
//		     fread(&t,sizeof(int),1,f);
//		     if((err = ferror(f)) != 0)
//		    	 {
//		    	 	 return err ;
//		    	 }
//
//	         *dbg_x = (double *)malloc(sizeof(double)*total_particles);
//	         debugPrintParticleCharacteristicArray(*dbg_x,total_particles,nt,"x",sort);
//
//	         *dbg_y = (double *)malloc(sizeof(double)*total_particles);
//	         debugPrintParticleCharacteristicArray(*dbg_y,total_particles,nt,"y",sort);
//
//	         *dbg_z = (double *)malloc(sizeof(double)*total_particles);
//	         debugPrintParticleCharacteristicArray(*dbg_z,total_particles,nt,"z",sort);
//
//	         *dbg_px = (double *)malloc(sizeof(double)*total_particles);
//	         debugPrintParticleCharacteristicArray(*dbg_px,total_particles,nt,"px",sort);
//
//	         *dbg_py = (double *)malloc(sizeof(double)*total_particles);
//	         debugPrintParticleCharacteristicArray(*dbg_py,total_particles,nt,"py",sort);
//
//	         *dbg_pz = (double *)malloc(sizeof(double)*total_particles);
//	         debugPrintParticleCharacteristicArray(*dbg_pz,total_particles,nt,"pz",sort);
//
//		 	readFortranBinaryArray(f,*dbg_x);
//		 	readFortranBinaryArray(f,*dbg_y);
//		 	readFortranBinaryArray(f,*dbg_z);
//		 	readFortranBinaryArray(f,*dbg_px);
//		 	readFortranBinaryArray(f,*dbg_py);
//		 	readFortranBinaryArray(f,*dbg_pz);
//
//		 	*qq_m = q_m;
//		 	*mm   = m;
//
//		 	if((err = ferror(f)) != 0)
//            {
//	   	 	    return err ;
//			}
//
//		 	return total_particles;
//      }

//	  virtual void readBinaryParticlesOneSort1(FILE *f,vector<Particle>& vp,
//			                                  particle_sorts sort,int nt)
//
//	  {
//		     double x,y,z,px,py,pz,q_m,m;
//		     int n = 0;
//		     Cell<Particle,dims> c0 = (*AllCells)[0];
//		     int pn_min,pn_ave,pn_max,pn_sum,err;
//
//
//		     if((err = ferror(f)) != 0) return;
//
//		     total_particles = readBinaryParticleArraysOneSort(f,&dbg_x,&dbg_y,&dbg_z,
//		    		                                             &dbg_px,&dbg_py,&dbg_pz,&q_m,&m,nt,
//		    		                                             sort);
//
//		     real_number_of_particle[(int)sort] = total_particles;
//
//		    err = ferror(f);
//		    for(int i = 0; i < total_particles;i++)
//		     {
//
//
//		    	  x   = dbg_x[i];
//		          y   = dbg_y[i];
//		          z   = dbg_z[i];
//  		          px   = dbg_px[i];
//		          py   = dbg_py[i];
//		          pz   = dbg_pz[i];
//
//
//			      Particle p;
//			      p.X.x   = x;
//			      p.X.y   = y;
//			      p.X.setZ(z);
//			      p.pu  = px;
//			      p.pv  = py;
//			      p.pw  = pz;
//			      p.m   = m;
//			      p.q_m = q_m;
//
//			      p.fortran_number = i+1;
//			      p.sort = sort;
//			      double3 d;
//			      d.x = x;
//			      d.y = y;
//			      d.z = z;
//
//			      n = c0.getPointCell(d);
//
//			      Cell<Particle,dims> & c = (*AllCells)[n];
//		    	  if((err = ferror(f)) != 0)
//		    	  {
//		    		  int qq = 0;
//		    	  }
//
//		    	  if(i == 3189003 && sort == PLASMA_ELECTRON)
//		    	  {
//		    		  int qq = 0;
//		    	  }
//
//		    	  printf("rank %d particle %d (%e,%e,%e) \n",getRank(),i+1,x,y,z);
//
//                  if(c.isPointInCell(p.GetX()) && (((i+1) % getSize()) == getRank()))
//                  {
//   			         if( (c.Insert(p) == true) )
//			         {
//#ifdef PARTICLE_PRINTS1000
//
//
//		        	          printf("rank %d particle %d (%e,%e,%e) is number %d in cell (%d,%d,%d)\n",
//		        	        		  getRank(),
//		        	    		 i+1,
//				    		       x,y,z,c.number_of_particles,c.cnum.x,c.cnum.y,c.cnum.z());
//
//
//#endif
//			         }
//			      }
//		     }// END total_particles LOOP
//
//		    err = ferror(f);
//		    free(dbg_x);
//			free(dbg_y);
//			free(dbg_z);
//		    free(dbg_px);
//			free(dbg_py);
//			free(dbg_pz);
//			err = ferror(f);
//
//             pn_min = 1000000000;
//             pn_max = 0;
//             pn_ave = 0;
//		     for(int n = 0;n < (*AllCells).size();n++)
//		     {
//		    	 Cell<Particle,dims> & c = (*AllCells)[n];
//
//		    	 pn_ave += c.number_of_particles;
//		    	 if(pn_min > c.number_of_particles) pn_min = c.number_of_particles;
//		    	 if(pn_max < c.number_of_particles) pn_max = c.number_of_particles;
//
//		     }
//		     err = ferror(f);
//		     pn_sum = pn_ave;
//		     pn_ave /= (*AllCells).size();
//
//		     printf("SORT m %15.5e q_m %15.5e %10d (sum %10d) particles in %8d cells: MIN %10d MAX %10d average %10d \n",
//		    		 m,            q_m,       total_particles,pn_sum,
//		    		 (*AllCells).size(),
//		    		 pn_min,pn_max,pn_ave);
//		     if((err = ferror(f)) != 0)
//		    		    	  {
//		    		    		 // int qq = 0;
//		    		    	  }
//
//		     err = ferror(f);
//	  }


	  virtual void readBinaryParticlesOneSort(FILE *f,vector<Particle>& vp,
			                                  particle_sorts sort,int nt,
			                                       		  double *dbg_x,
			                                      		  double *dbg_y,
			                                      		  double *dbg_z,
			                                      		  double *dbg_px,
			                                      		  double *dbg_py,
			                                      		  double *dbg_pz,
			                                      		  double q_m,
			                                      		  double m,
			                                      		  double *masses,
			                                      		  int total_particles,
			                                      		  int add_flag
			                                  )

	  {
		     double x,y,z,px,py,pz;//,q_m,m;
		     int n = 0;
		     Cell<Particle,dims> c0 = (*AllCells)[0];
		     int pn_min,pn_ave,pn_max,pn_sum,err;
		     static int first = 1;


//		     if((err = ferror(f)) != 0) return;

             if(first == 1)
             {
            	 real_number_of_particle[0] = 0;
            	 real_number_of_particle[1] = 0;
            	 real_number_of_particle[2] = 0;
            	 first = 0;
             }

		     real_number_of_particle[(int)sort] += total_particles;

//		    err = ferror(f);
		    for(int i = 0; i < total_particles;i++)
		     {

		          
		    	  x   = dbg_x[i];
		          y   = dbg_y[i];
		          z   = dbg_z[i];
  		          px   = dbg_px[i];
		          py   = dbg_py[i];
		          pz   = dbg_pz[i];


			      Particle p;
			      p.X.x   = x;
			      p.X.y   = y;
			      p.X.setZ(z);
			      p.pu  = px;
			      p.pv  = py;
			      p.pw  = pz;
			      if(beam_plasma == 1)
			      {
			         p.m   = m;
			      }
			      else
			      {
  			    	  p.m  = masses[i];
			      }
			      p.q_m = q_m;

			      p.fortran_number = i+1;
			      p.sort = sort;
			      double3 d;
			      d.x = x;
			      d.y = y;
			      d.z = z;

			      n = c0.getPointCell(d);

			      Cell<Particle,dims> & c = (*AllCells)[n];
//		    	  if((err = ferror(f)) != 0)
//		    	  {
//		    		  int qq = 0;
//		    	  }

		    	  if(i == 3189003 && sort == PLASMA_ELECTRON)
		    	  {
		    		  int qq = 0;
		    	  }

                  if(c.isPointInCell(p.GetX()) && (((i+1) % getSize()) == getRank()))
                  {
   			         if( (c.Insert(p) == true) )
			         {
   			        	 if(add_flag == 1)
   			        	 {
#ifdef PARTICLE_PRINTS1000
   			        		 printf("rank %d particle %d (%e,%e,%e) is number %d in cell (%d,%d,%d)\n",getRank(),
   			        				        	    		 i+1,
   			        						    		       x,y,z,c.number_of_particles,c.cnum.x,c.cnum.y,c.cnum.z());
#endif
   			        	 }
#ifdef PARTICLE_PRINTS1000
//		                  if((i+1) == 0 )
//		                  {
		        	          printf("rank %d particle %d (%e,%e,%e) is number %d in cell (%d,%d,%d)\n",getRank(),
		        	    		 i+1,
				    		       x,y,z,c.number_of_particles,c.cnum.x,c.cnum.y,c.cnum.z());
//		                  }

#endif
			         }
			      }
		     }// END total_particles LOOP

//		    err = ferror(f);
		    free(dbg_x);
			free(dbg_y);
			free(dbg_z);
		    free(dbg_px);
			free(dbg_py);
			free(dbg_pz);
//			err = ferror(f);

             pn_min = 1000000000;
             pn_max = 0;
             pn_ave = 0;
		     for(int n = 0;n < (*AllCells).size();n++)
		     {
		    	 Cell<Particle,dims> & c = (*AllCells)[n];

		    	 pn_ave += c.number_of_particles;
		    	 if(pn_min > c.number_of_particles) pn_min = c.number_of_particles;
		    	 if(pn_max < c.number_of_particles) pn_max = c.number_of_particles;

		     }
//		     err = ferror(f);
		     pn_sum = pn_ave;
		     pn_ave /= (*AllCells).size();

		     printf("SORT m %15.5e q_m %15.5e %10d (sum %10d) particles in %8lu cells: MIN %10d MAX %10d average %10d \n",
		    		 m,      
	                        q_m,       
	                         total_particles,
	                          pn_sum,
		    		 (*AllCells).size(),
		    		 pn_min,
	                         pn_max,
	                         pn_ave);
//		     if((err = ferror(f)) != 0)
//		    		    	  {
//		    		    		 // int qq = 0;
//		    		    	  }
//
//		     err = ferror(f);
	  }

//	  void InitBinaryParticlesArrays(char *fn,int nt,
//	     double **dbg_ion_x, double **dbg_ion_y, double **dbg_ion_z,
//	     double **dbg_ion_px,double **dbg_ion_py,double **dbg_ion_pz,double *ion_q_m,double *ion_m,
//	     int *total_ions,
//	     double **dbg_el_x, double **dbg_el_y, double **dbg_el_z,
//	     double **dbg_el_px,double **dbg_el_py,double **dbg_el_pz,double *el_q_m,double *el_m,
//	     int *total_electrons,
//	     double **dbg_beam_x, double **dbg_beam_y, double **dbg_beam_z,
//	     double **dbg_beam_px,double **dbg_beam_py,double **dbg_beam_pz,
//	     double *beam_q_m,double *beam_m,
//	     int *total_beam_electrons)
//	  {
//	     FILE *f;
//	     double *buf;
//
//
//
//	     buf = (double *)malloc(sizeof(double)*mesh.size2());
//
//	     if((f = fopen(fn,"rb")) == NULL) return;
//	     struct sysinfo info;
//
//	     sysinfo(&info);
//	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
//	     readFortranBinaryArray(f,buf);
//	     sysinfo(&info);
//	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
//	     readFortranBinaryArray(f,buf);
//	     sysinfo(&info);
//	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
//
//	     readFortranBinaryArray(f,buf);
//	     sysinfo(&info);
//	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
//
//	     readFortranBinaryArray(f,buf);
//	     sysinfo(&info);
//	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
//
//	     readFortranBinaryArray(f,buf);
//	     readFortranBinaryArray(f,buf);
//
//	     readFortranBinaryArray(f,buf);
//	     readFortranBinaryArray(f,buf);
//	     readFortranBinaryArray(f,buf);
//
//	     readFortranBinaryArray(f,buf);
//	     readFortranBinaryArray(f,buf);
//	     readFortranBinaryArray(f,buf);
//	     int err;
//	     err = ferror(f);
//
//	     *total_ions = readBinaryParticleArraysOneSort(f,     dbg_ion_x, dbg_ion_y, dbg_ion_z,
//	    		                                             dbg_ion_px,dbg_ion_py,dbg_ion_pz,
//	    		                                             ion_q_m,ion_m,nt,ION);
//
////	     readBinaryParticlesOneSort(f,vp,ION,nt,dbg_ion_x,dbg_ion_y,dbg_ion_z,
////                                               dbg_ion_px,dbg_ion_py,dbg_ion_pz,
////                                               ion_q_m,ion_m,total_ions);
//
//	     sysinfo(&info);
//	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
//	     err = ferror(f);
//
////	     readBinaryParticlesOneSort1(f,vp,PLASMA_ELECTRON,nt);
//	     *total_electrons = readBinaryParticleArraysOneSort(f,     dbg_el_x, dbg_el_y, dbg_el_z,
//	    		                                             dbg_el_px,dbg_el_py,dbg_el_pz,
//	    		                                             el_q_m,el_m,nt,PLASMA_ELECTRON);
//
////	     readBinaryParticlesOneSort(f,vp,PLASMA_ELECTRON,nt,dbg_el_x,dbg_el_y,dbg_el_z,
////                                               dbg_el_px,dbg_el_py,dbg_el_pz,
////                                               el_q_m,el_m,total_electrons);
//
//
//	     err = ferror(f);
//	     sysinfo(&info);
//	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
//         err = ferror(f);
////	     readBinaryParticlesOneSort1(f,vp,BEAM_ELECTRON,nt);
//	     *total_beam_electrons = readBinaryParticleArraysOneSort(f,dbg_beam_x,dbg_beam_y,dbg_beam_z,
//	    		                                             dbg_beam_px,dbg_beam_py,dbg_beam_pz,
//	    		                                             beam_q_m,beam_m,nt,BEAM_ELECTRON);
//
////	     readBinaryParticlesOneSort(f,vp,BEAM_ELECTRON,nt,dbg_beam_x,dbg_beam_y,dbg_beam_z,
////                                               dbg_beam_px,dbg_beam_py,dbg_beam_pz,
////                                               beam_q_m,beam_m,total_beam_electrons);
//
//	     fclose(f);
//
//	     magf = 1;
//	  }

	  virtual void InitBinaryParticles(vector<Particle>& vp,int nt,
			  double *dbg_ion_x, double *dbg_ion_y, double *dbg_ion_z,
			  double *dbg_ion_px,double *dbg_ion_py,double *dbg_ion_pz,
			  double ion_q_m,double ion_m,
			  int total_ions,
			  double *dbg_el_x, double *dbg_el_y, double *dbg_el_z,
			  double *dbg_el_px,double *dbg_el_py,double *dbg_el_pz,
			  double el_q_m,double el_m,
			  int total_electrons,
			  double *dbg_beam_x, double *dbg_beam_y, double *dbg_beam_z,
			  double *dbg_beam_px,double *dbg_beam_py,double *dbg_beam_pz,
			  double beam_q_m,double beam_m,
			  int total_beam_electrons
	  )
	  {
	     FILE *f;
         double *masses;
	     readBinaryParticlesOneSort(f,vp,ION,nt,dbg_ion_x,dbg_ion_y,dbg_ion_z,
                                               dbg_ion_px,dbg_ion_py,dbg_ion_pz,
                                               ion_q_m,ion_m,
                                               masses,total_ions,0);

//	     sysinfo(&info);
//	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
//	     err = ferror(f);

//	     readBinaryParticlesOneSort1(f,vp,PLASMA_ELECTRON,nt);
//	     total_electrons = readBinaryParticleArraysOneSort(f,     &dbg_el_x, &dbg_el_y, &dbg_el_z,
//	    		                                             &dbg_el_px,&dbg_el_py,&dbg_el_pz,
//	    		                                             &el_q_m,&el_m,nt,PLASMA_ELECTRON);

	     readBinaryParticlesOneSort(f,vp,PLASMA_ELECTRON,nt,dbg_el_x,dbg_el_y,dbg_el_z,
                                               dbg_el_px,dbg_el_py,dbg_el_pz,
                                               el_q_m,el_m,masses,
                                               total_electrons,0);


//	     err = ferror(f);
//	     sysinfo(&info);
//	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
//         err = ferror(f);
//	     readBinaryParticlesOneSort1(f,vp,BEAM_ELECTRON,nt);
//	     total_beam_electrons = readBinaryParticleArraysOneSort(f,     &dbg_beam_x, &dbg_beam_y, &dbg_beam_z,
//	    		                                             &dbg_beam_px,&dbg_beam_py,&dbg_beam_pz,
//	    		                                             &beam_q_m,&beam_m,nt,BEAM_ELECTRON);

	     readBinaryParticlesOneSort(f,vp,BEAM_ELECTRON,nt,dbg_beam_x,dbg_beam_y,dbg_beam_z,
                                               dbg_beam_px,dbg_beam_py,dbg_beam_pz,
                                               beam_q_m,beam_m,masses,
                                               total_beam_electrons,0);

//	     fclose(f);

	     magf = 1;
	  }

	  virtual void InitBinaryParticlesTotal(vector<Particle>& vp,int nt,
			  double *dbg_ion_x, double *dbg_ion_y, double *dbg_ion_z,
			  double *dbg_ion_px,double *dbg_ion_py,double *dbg_ion_pz,
			  double ion_q_m,double *ion_m,
			  int total_ions
	  )
	  {
	     FILE *f;


	     double m = 1.0; // TODO: FIX IT
	     readBinaryParticlesOneSort(f,vp,ION,nt,dbg_ion_x,dbg_ion_y,dbg_ion_z,
                                               dbg_ion_px,dbg_ion_py,dbg_ion_pz,
                                               ion_q_m,m,ion_m,total_ions,0);



//	     fclose(f);

	     magf = 1;
	  }


	  virtual void InitElectronParticles(){}
//	  virtual void InitIonParticles(int n_per_cell1,double q_m,vector<Particle> &vecp)
//	  {
//	     int total_ions = mesh.x*mesh.y*mesh.z()*n_per_cell;
//	     Particle *p;
//	     double x,y,z;
//
//	     for(int j = 0;j < total_ions;j++)
//	     {
//		z = Lz * rnd_uniform();
//		y = meh * Ly + Ly * rnd_uniform();
//		x = Lx * rnd_uniform();
//
//		p = new Particle(x,y,z,0.0,0.0,0.0,ni,q_m);
//
//	#ifdef DEBUG_PLASMA
//	#endif
//
//		vecp.push_back(*p);
//	     }
//	  }

	  virtual void InitBeamParticles(int n_per_cell1){}
	  void Distribute(vector<Particle> &vecp)
	  {
	     Cell<Particle,dims> c0 = (*AllCells)[0],c111;
	     int    n;//,i;
	     int  vec_size = vecp.size();

	     for(int j = 0;j < vecp.size();j++)
	     {
		 Particle p = vecp[j];
		 double3 d;
		 d.x = p.X.x;
		 d.y = p.X.y;
		 d.z = p.X.z();

		 n = c0.getPointCell(d);

		 Cell<Particle,dims> & c = (*AllCells)[n];;

		 if(c.Insert(p) == true)
		 {
	#ifdef PARTICLE_PRINTS1000
         if((vec_size-vecp.size())%1000 == 0 )	printf("particle %d (%e,%e,%e) is number %d in cell (%d,%d,%d)\n",vec_size-vecp.size(),
		    		p.x,p.y,p.z,c.number_of_particles,c.i,c.l,c.k);
         if((vec_size-vecp.size()) == 10000) exit(0);
	#endif
		    vecp.erase(vecp.begin()+j);
		    j--;
		 }
	     }
	     int pn_min = 1000000,pn_max = 0,pn_ave = 0;
	     for(int n = 0;n < (*AllCells).size();n++)
	     {
	    	 Cell<Particle,dims> & c = (*AllCells)[n];

	    	 pn_ave += c.number_of_particles;
	    	 if(pn_min > c.number_of_particles) pn_min = c.number_of_particles;
	    	 if(pn_max < c.number_of_particles) pn_max = c.number_of_particles;

	     }
	     pn_ave /= (*AllCells).size();

	     printf("%10d particles in %8d cells: MIN %5d MAX %5d average %5d \n",vec_size,(*AllCells).size(),
	    		                                                              pn_min,pn_max,pn_ave);
	  }



	void virtual emeIterate(int i_s,int i_f,int l_s,int l_f,int k_s,int k_f,
			double *E,double *H1, double *H2,
			double *J,double c1,double c2, double tau,
			int dx1,int dy1,int dz1,int dx2,int dy2,int dz2)
	{
		Cell<Particle,dims>  c0 = (*AllCells)[0],*c;

		c = &c0;

		for(int i = i_s;i <= i_f;i++)
		{
			  for(int l = l_s;l <= i_f;l++)
			  {
			      for(int k = k_s;k <= k_f;k++)
			      {
			    	  emeElement(c,i,l,k,E,H1,H2,
			    	  		J,c1,c2,tau,
			    	  		dx1,dy1,dz1,dx2,dy2,dz2);
			      }
			  }
		}
	}


	int MagneticFieldTrace(Cell<Particle,dims> &c,string lname,int nt,double *Q,double *H,double *E1,double *E2,int i_end,int l_end,int k_end,double c1,double c2,int dir)
	{
	      int dx1,dy1,dz1,dy2,dx2,dz2;
	      double e1_n1,e1_n,e2_n2,e2_n,t;

	#ifdef DEBUG_PLASMA_FIELDS
	{    char logname[100];
	     double *Hres,*ldQ,*ldH,*ldE1,*ldE2;

	     Hres = (double *)malloc(mesh.size2()*sizeof(double));
	     ldQ  = (double *)malloc(mesh.size2()*sizeof(double));
	     ldE1 = (double *)malloc(mesh.size2()*sizeof(double));
	     ldE2 = (double *)malloc(mesh.size2()*sizeof(double));
	     ldH  = (double *)malloc(mesh.size2()*sizeof(double));

	     sprintf(logname,"%s%03d.dat",lname.c_str(),nt);

	     read3DarrayLog(logname, ldQ,40,0);
	     read3DarrayLog(logname, ldE1,40,2);
	     CheckArray(E1,ldE1);
	     read3DarrayLog(logname, ldE2,40,4);
	     CheckArray(E2,ldE2);

	     read3DarrayLog(logname, ldH,40,5);
	     CheckArray(H,ldH);
	     n0 = c.getGlobalCellNumber(11,10,10);
	     printf("%e \n",H[n0]);
	     read3DarrayLog(logname, Hres,40,6);
	     free(Hres);
	     free(ldQ);
	     free(ldH);
	     free(ldE1);
	     free(ldE2);
	}
	#endif

	      dx1 = (dir == 0)*0 + (dir == 1)*1 + (dir == 2)*0;
	      dy1 = (dir == 0)*0 + (dir == 1)*0 + (dir == 2)*1;
	      dz1 = (dir == 0)*1 + (dir == 1)*0 + (dir == 2)*0;

	      dx2 = (dir == 0)*0 + (dir == 1)*0 + (dir == 2)*1;
	      dy2 = (dir == 0)*1 + (dir == 1)*0 + (dir == 2)*0;
	      dz2 = (dir == 0)*0 + (dir == 1)*1 + (dir == 2)*0;

     if(CPU_field == 0)
     {
//   		dim3 dimGrid(i_end+1,l_end+1,k_end+1),dimBlock(1,1,1);
   		int err;


   		params.i_s = 0;
   		params.k_s = 0;
   		params.l_s = 0;
   		params.H = H;
   		params.Q = Q;
   		params.E1 = E1;
   		params.E2 = E2;
   		params.c1 = c1;
   		params.c2 = c2;
   		params.dx1 = dx1;
   		params.dy1 = dy1;
   		params.dz1 = dz1;
   		params.dx2 = dx2;
   		params.dy2 = dy2;
   		params.dz2 = dz2;

   		MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

   		Kernel_Launcher(d_CellArray,d_params,
   				   	    i_end+1,l_end+1,k_end+1,
   				   			    		1,1,1,
   				   			    		16000,
   				   			    		h_emh1,
   				   			    		"emh1");

//	    GPU_emh1<<<dimGrid,dimBlock>>>(d_CellArray,0,0,0,Q,H,E1,E2,c1,c2,
//	    		dx1,dy1,dz1,dx2,dy2,dz2);
	    err = getLastError();
	    ThreadSynchronize();
//	    TEST_ERROR(err);
     }
     else
     {

	      for(int k = 0;k <= k_end;k++)
	      {
	   	    for(int l = 0;l <= l_end;l++)
		    {
		      for(int i = 0;i <= i_end;i++)
		      {
			  int n  = c.getGlobalCellNumber(i,l,k);
			  int n1 = c.getGlobalCellNumber(i+dx1,l+dy1,k+dz1);
			  int n2 = c.getGlobalCellNumber(i+dx2,l+dy2,k+dz2);

			  e1_n1 = E1[n1];
			  e1_n  = E1[n];
			  e2_n2 = E2[n2];
			  e2_n  = E2[n];

			  t  = 0.5*(c1*(e1_n1 - e1_n)- c2*(e2_n2 - e2_n));
	#ifdef FIELDS_DEBUG_OUTPUT
			  std::cout << i << l << k << t << " " << Q[n] << " diff " << fabs(t - Q[n]) << std::endl;

			  std::cout << "E1 " <<  E1[n] << " E2 " << E2[n] << " H  " << H[n] << std::endl;

			  if((fabs(t - Q[n]) > TOLERANCE) || ((i == 0) && (l == 0) && (k == 10)))
			  {
			     printf("i %3d l %3d k %3d %15.5e dbg %15.5e\n",i+1,l+1,k+1,t,Q[n]);
			  }
	#endif
			  Q[n] = t;
#ifdef FIELDS_DEBUG_OUTPUT
			  printf("%d %d %d %e %e \n",i,l,k,H[n],Q[n]);
#endif
			  H[n] += Q[n];
		      }
		    }
	      }
     }
     return 0;
	  }

	int AddConstantMagneticField(Cell<Particle,dims> &c,double *H,int i_end,int l_end,int k_end,double B0,int nt)
	{
		int err;
//		if(CPU_field == 0)
//		{
//		   		dim3 dimGrid(i_end+1,l_end+1,k_end+1),dimBlock(1,1,1);
		   		// int i_s,int l_s,int k_s,
//		   		params.Q = Q;
//        		printGPUArray	(H,502,nt,"Bx");
		   		params.H = H;
		   		params.i_s = 0;
		   		params.l_s = 0;
		   		params.k_s = 0;
		   		params.B0  = B0;

		   		MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

		   	    Kernel_Launcher(d_CellArray,d_params,
		   	    		        i_end+1,l_end+1,k_end,
		   			    		1,1,1,
		   			    		16000,
		   			    		h_add,
		   			    		"addB");
//		   	    printGPUArray	(H,509,nt,"Bx");

//			    GPU_emh2<<<dimGrid,dimBlock>>>(d_CellArray,0,0,0,Q,H);
			    err = getLastError();
			    ThreadSynchronize();
//			    TEST_ERROR(err);
//		}

	}

	int SimpleMagneticFieldTrace(Cell<Particle,dims> &c,double *Q,double *H,int i_end,int l_end,int k_end)
	{
		int err;
		if(CPU_field == 0)
		{
//		   		dim3 dimGrid(i_end+1,l_end+1,k_end+1),dimBlock(1,1,1);
		   		// int i_s,int l_s,int k_s,
		   		params.Q = Q;
		   		params.H = H;
		   		params.i_s = 0;
		   		params.l_s = 0;
		   		params.k_s = 0;

		   		MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

		   	    Kernel_Launcher(d_CellArray,d_params,
		   	    		        i_end+1,l_end+1,k_end+1,
		   			    		1,1,1,
		   			    		16000,
		   			    		h_emh2,
		   			    		"emh2");

//			    GPU_emh2<<<dimGrid,dimBlock>>>(d_CellArray,0,0,0,Q,H);
			    err = getLastError();
			    ThreadSynchronize();
//			    TEST_ERROR(err);
		}
		else
		{
	      for(int i = 0;i <= i_end;i++)
	      {
		  for(int l = 0;l <= l_end;l++)
		  {
		      for(int k = 0;k <= k_end;k++)
		      {
			      int n = c.getGlobalCellNumber(i,l,k);
			      H[n] += Q[n];
		      }
		  }
	      }
		}

	      return 0;
	  }
	  int PeriodicBoundaries(double *E,int dir,int start1,int end1,int start2,int end2,int N)
	  {
	      Cell<Particle,dims>  c = (*AllCells)[0];
	      int err;

	      if(CPU_field == 0)
	      {
//	    		dim3 dimGrid(end1-start1+1,1,end2-start2+1),dimBlock(1,1,1);

	    		params.E    = E;
	    		params.dir  = dir;
	    		params.i_s  = start1;
	    		params.k_s  = start2;
	    		params.to   = 0;
	    		params.from = N;

	    		MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

	    		Kernel_Launcher(d_CellArray,d_params,
	    				end1-start1+1,1,end2-start2+1,
	    		    		             1,1,1,
	    		    		             16000,
	    		                         h_GPU_periodic_SingleNode,
	    		                         "periodic");

//	    	    GPU_periodic<<<dimGrid,dimBlock>>>(d_CellArray,start1,start2,E,dir,0,N);
	    	    err = getLastError();
	    	    ThreadSynchronize();
//	    	    TEST_ERROR(err);

	    	    params.to   = N+1;
	    	    params.from = 1;

	    		MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

	    		Kernel_Launcher(d_CellArray,d_params,
	    				end1-start1+1,1,end2-start2+1,
	    		    		             1,1,1,
	    		    		             16000,
	    		                         h_GPU_periodic_SingleNode,
	    		                         "periodic");



//	    	    GPU_periodic<<<dimGrid,dimBlock>>>(d_CellArray,start1,start2,E,dir,N+1,1);
	    	    err = getLastError();
	    	    ThreadSynchronize();
//	    	    TEST_ERROR(err);

	      }
	      else
	      {

	         for(int k = start2;k <= end2;k++)
	         {
		        for(int i = start1;i <= end1;i++)
		        {
			       periodicElement(&c,i,k,E,dir,0,N);
		        }
	         }
	         for(int k = start2;k <= end2;k++)
	         {
	            for(int i = start1;i <= end1;i++)
	      	    {
	        	    periodicElement(&c,i,k,E,dir,N+1,1);
		        }
	         }
	      }
	      return 0;
	}

int SinglePeriodicBoundary(double *E,int dir,int start1,int end1,int start2,int end2,int N)
{
    Cell<Particle,dims>  c = (*AllCells)[0];

    if(CPU_field == 0)
    {
    	int err;
//    	dim3 dimGrid(end1-start1+1,1,end2-start2+1),dimBlock(1,1,1);

		params.E    = E;
		params.dir  = dir;
		params.i_s  = start1;
		params.k_s  = start2;
		params.to   = N+1;
		params.from = 1;

		MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

		Kernel_Launcher(d_CellArray,d_params,
				end1-start1+1,1,end2-start2+1,
		    		             1,1,1,
		    		             16000,
		                         h_GPU_periodic_SingleNode,
		                         "periodic");

//   	    GPU_periodic<<<dimGrid,dimBlock>>>(d_CellArray,start1,start2,E,dir,N+1,1);
   	    err = getLastError();
   	    ThreadSynchronize();
//   	    TEST_ERROR(err);

    }
    else
    {
       for(int k = start2;k <= end2;k++)
       {
	  	  for(int i = start1;i <= end1;i++)
	  	  {
	  		  int3 i0,i1;

              int n   = c.getGlobalBoundaryCellNumber(i,k,dir,N+1);
              int n1  = c.getGlobalBoundaryCellNumber(i,k,dir,1);
	  		  E[n]    = E[n1];
	  		  i0      = c.getCellTripletNumber(n);
	  		  i1      = c.getCellTripletNumber(n1);
	  		  std::cout << "ex1 "<< i0.x+1 << " "<< i0.y+1 << " " << i0.z+1  <<" " <<  i1.x+1 << " " << i1.y+1 << " " << i1.z+1  << " " << E[n]  << " " << E[n1] << std::endl;
	   	  }
        }
    }
    return 0;
}



      void getWrongCurrentCellList(int num,int nt)
      {
    	   FILE *f;
    	   char fn_copy[100];
    	   Cell<Particle,dims> c = (*AllCells)[0];
    	   int wrong_flag[mesh.size2()],*d_wrong_flag;
    	   double_pointer wrong_attributes[mesh.size2()],*d_wrong_attributes;

    	   readControlPoint(&f,fn_copy,num,nt,1,0,dbgEx,dbgEy,dbgEz,dbgHx,dbgHy,dbgHz,dbgJx,dbgJy,dbgJz,dbg_Qx,dbg_Qy,dbg_Qz,
    	                     dbg_x,dbg_y,dbg_z,dbg_px,dbg_py,dbg_pz);

    	   for(int n = 0;n < mesh.size2();n++)
    	       	   {
    		           wrong_flag[n] = 0;
    	       	   }
    	   MemoryAllocate((void**)&d_wrong_flag,sizeof(int)*mesh.size2());
    	   MemoryAllocate((void**)&d_wrong_attributes,sizeof(double_pointer)*mesh.size2());

    	   static double *t;
    	   static int first = 1;

    	   if(first == 1)
    	   {
    	  	 t = (double *)malloc(sizeof(double)*mesh.size2());
    	  	 first = 0;
    	   }
    	   int err;
    	   err = MemoryCopy(t,d_Jx,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
    	   if(err != 0)
    	   {
    	     printf("getWrongCurrentCellList err %d %s \n",err,getErrorString(err));
    	  	 exit(0);
    	   }

    	   double diff = 0.0;
    	   jx_wrong_points_number = 0;
    	   for(int n = 0;n < mesh.size2();n++)
    	   {
   	           if(fabs(dbgJx[n] - t[n]) > TOLERANCE)
    	   	   {
      		       int3 i = c.getCellTripletNumber(n);
                   jx_wrong_points_number++;
                   wrong_flag[n] = 1;
        	   }
    	   }
    	   jx_wrong_points = (int3 *)malloc(jx_wrong_points_number*sizeof(int3));

    	   int num_cell = 0;
    	   for(int n = 0;n < mesh.size2();n++)
    	   {
    	       if(fabs(dbgJx[n] - t[n]) > TOLERANCE)
    	       {
    	          int3 i = c.getCellTripletNumber(n);
    	          jx_wrong_points[num_cell++] = i;
    	       }
    	   }
    	   MemoryAllocate((void**)&(d_jx_wrong_points),jx_wrong_points_number*sizeof(int3));

    	   MemoryCopy(d_wrong_flag,wrong_flag,sizeof(int)*mesh.size2(),HOST_TO_DEVICE);

    	   MemoryCopy(d_wrong_attributes,wrong_attributes,
    			                              sizeof(double_pointer)*mesh.size2(),HOST_TO_DEVICE);

//           copy_pointers<<<mesh.size2(),1>>>(d_CellArray,d_wrong_flag,d_wrong_attributes);
           err = getLastError();
           ThreadSynchronize();
//           TEST_ERROR(err);

      }

      void WrongCurrentCell_AttributeMalloc(int num,int nt)
      {

      }

	  double checkFirstHalfstepFields(int nt)
	  {
		  return 1.0;//(t_ex+t_ey+t_ez+t_hx+t_hy+t_hz)/6.0;
	  }

	  double checkFirstHalfstep_emh2_GPUMagneticFields(int nt)
	  	  {
	  		  double /*t = 0.0,*/*dbg/*,t_ex,t_ey,t_ez*/,t_hx,t_hy,t_hz;

	  		  dbg = (double *)malloc(sizeof(double)*mesh.size2());

	  		  return 1.0;
	  	  }


	  double checkFirstHalfstep_emh1_GPUMagneticFields(int nt)
	  {
		  double /*t = 0.0,*/*dbg,t_ex,t_ey,t_ez,t_hx,t_hy,t_hz;

		  return 1.0;//(t_ex+t_ey+t_ez+t_hx+t_hy+t_hz)/6.0;
	  }

	  double checkFirstHalfstep_emh1_MagneticFields(int nt,double *Qx,double *Qy,double *Qz,
			                                               double *Hx,double *Hy,double *Hz)
	  {
		  double *dbg,t_ex,t_ey,t_ez,t_hx,t_hy,t_hz;

		  return 1.0;//(t_ex+t_ey+t_ez+t_hx+t_hy+t_hz)/6.0;
	  }


	  double checkFirstHalfstepElectricFields(int nt)
	  {
		  double *dbg,t_ex,t_ey,t_ez;//,t_hx,t_hy,t_hz;

		  return 1.0;//(t_ex+t_ey+t_ez)/3.0;
	  }

	  double checkFirstHalfstepGPUElectricFields(int nt)
	  {
		  return 1.0;//(t_ex+t_ey+t_ez)/3.0;
	  }

	  double checkSecondHalfstepFields(int nt)
	    {
	  	  return 1.0;//(t_ex+t_ey+t_ez)/3.0;
	    }
	  double checkGPUSecondHalfstepFields(int nt)
	  	    {
	  	  	  return 1.0;//(t_ex+t_ey+t_ez)/3.0;
	  	    }


	  double checkCurrents(int nt)
	    {
	  	  double  *dbg,/*t_ex,t_ey,t_ez,*/t_hx,t_hy,t_hz;
	  	  static int first = 1;

	  	  if(first == 1)
	  	  {
	  	     dbg = (double *)malloc(sizeof(double)*mesh.size2());
	  	     first = 0;
	  	  }

	  	  return 1.0;//(t_hx+t_hy+t_hz)/3.0;
	    }

	  void SetPeriodicCurrents(int nt)
	  {
          double *dbg = (double *)malloc(mesh.size2()*sizeof(double));
          int err;

          checkGPUArray(Jx,d_Jx);
	      PeriodicCurrentBoundaries(Jx,dbgJx,0,0, 0,mesh.y+1, 0,mesh.z()+1);

	      dim3 dimGridX(mesh.y+2,1,mesh.dimz2()),
	    	   dimGridY(mesh.x+2,1,mesh.dimz2()),dimGridZ(mesh.x+2,1,mesh.y+2),
	    	   dimBlock(1,1,1);

	      int N = getBoundaryLimit(0);

	      //double *E,int dirE, int dir,
          //int i_s,int k_s,int N
	      params.E = d_Jx;
	      params.dirE = 0;
	      params.dir  = 0;
	      params.i_s  = 0;
	      params.k_s  = 0;
	      params.N    = mesh.x+2;
	      MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

	  	  Kernel_Launcher(d_CellArray,d_params,dimGridX.x,dimGridX.y,dimGridX.z,
	  	                dimBlock.x,dimBlock.y,dimBlock.z,16000,h_CurrentPeriodic,
	  	               "CurrentPeriodic");

//	      GPU_CurrentPeriodic<<<dimGridX,dimBlock>>>(d_CellArray,d_Jx,0,0,0,0,mesh.x+2);
	      err = getLastError();
	      ThreadSynchronize();
//	      TEST_ERROR(err);


	     checkGPUArray(Jx,d_Jx);
	     PeriodicCurrentBoundaries(Jx,dbgJx,0,1,0,mesh.x+1,0,mesh.z()+1);
	      params.E = d_Jx;
	      params.dirE = 0;
	      params.dir  = 1;
	      params.i_s  = 0;
	      params.k_s  = 0;
	      params.N    = mesh.y+2;
	      MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

	  	  Kernel_Launcher(d_CellArray,d_params,dimGridY.x,dimGridY.y,dimGridY.z,
	  	                dimBlock.x,dimBlock.y,dimBlock.z,16000,h_CurrentPeriodic,
	  	               "CurrentPeriodic");

//	     GPU_CurrentPeriodic<<<dimGridY,dimBlock>>>(d_CellArray,d_Jx,0,1,0,0,mesh.y+2);
	     err = getLastError();
	     ThreadSynchronize();
//	     TEST_ERROR(err);
	     checkGPUArray(Jx, d_Jx);
	     PeriodicCurrentBoundaries(Jx,dbgJx,0,2,0,mesh.x+1,0,mesh.y+1);

	     params.E = d_Jx;
	     params.dirE = 0;
	     params.dir  = 2;
	     params.i_s  = 0;
	     params.k_s  = 0;
	     params.N    = mesh.dimz2();
	     MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

	     Kernel_Launcher(d_CellArray,d_params,dimGridZ.x,dimGridZ.y,dimGridZ.z,
	                     dimBlock.x,dimBlock.y,dimBlock.z,16000,h_CurrentPeriodic,
	                    "CurrentPeriodic");

//	     GPU_CurrentPeriodic<<<dimGridZ,dimBlock>>>(d_CellArray,d_Jx,0,2,0,0,mesh.dimz2());
	     err = getLastError();
	     ThreadSynchronize();
//	     TEST_ERROR(err);
	     	          checkGPUArray(Jx, d_Jx);

	     checkGPUArray(Jy, d_Jy);
	     PeriodicCurrentBoundaries(Jy,dbgJy,1,0,0,mesh.y+1,0,mesh.z()+1);

	     params.E    = d_Jy;
	     params.dirE = 1;
	     params.dir  = 0;
	     params.i_s  = 0;
	     params.k_s  = 0;
	     params.N    = mesh.x+2;
	     MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

	     Kernel_Launcher(d_CellArray,d_params,dimGridX.x,dimGridX.y,dimGridX.z,
	                     dimBlock.x,dimBlock.y,dimBlock.z,16000,h_CurrentPeriodic,
	                     "CurrentPeriodic");

//	     GPU_CurrentPeriodic<<<dimGridX,dimBlock>>>(d_CellArray,d_Jy,1,0,0,0,mesh.x+2);
	     err = getLastError();
	     ThreadSynchronize();
//	     TEST_ERROR(err);
	     checkGPUArray(Jy, d_Jy);
	     PeriodicCurrentBoundaries(Jy,dbgJy,1,1,0,mesh.x+1,0,mesh.z()+1);

	     params.E    = d_Jy;
	     params.dirE = 1;
	     params.dir  = 1;
	     params.i_s  = 0;
	     params.k_s  = 0;
	     params.N    = mesh.y+2;
	     MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

	     Kernel_Launcher(d_CellArray,d_params,dimGridY.x,dimGridY.y,dimGridY.z,
	                     dimBlock.x,dimBlock.y,dimBlock.z,16000,h_CurrentPeriodic,
	                     "CurrentPeriodic");

//	     GPU_CurrentPeriodic<<<dimGridY,dimBlock>>>(d_CellArray,d_Jy,1,1,0,0,mesh.y+2);
	     err = getLastError();
	     ThreadSynchronize();
//	     TEST_ERROR(err);
	     checkGPUArray(Jy, d_Jy);
	     PeriodicCurrentBoundaries(Jy,dbgJy,1,2,0,mesh.x+1,0,mesh.y+1);

	     params.E    = d_Jy;
  	     params.dirE = 1;
   	     params.dir  = 2;
  	     params.i_s  = 0;
   	     params.k_s  = 0;
   	     params.N    = mesh.dimz2();
   	     MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);
   	     Kernel_Launcher(d_CellArray,d_params,dimGridZ.x,dimGridZ.y,dimGridZ.z,
	     	                     dimBlock.x,dimBlock.y,dimBlock.z,16000,h_CurrentPeriodic,
	     	                     "CurrentPeriodic");

//	     GPU_CurrentPeriodic<<<dimGridZ,dimBlock>>>(d_CellArray,d_Jy,1,2,0,0,mesh.dimz2());
	     err = getLastError();
	     ThreadSynchronize();
//	     TEST_ERROR(err);
	     checkGPUArray(Jy, d_Jy);
	     checkGPUArray(Jz, d_Jz);
	     PeriodicCurrentBoundaries(Jz,dbgJz,2,0,0,mesh.y+1,0,mesh.z()+1);

	     params.E    = d_Jz;
   	     params.dirE = 2;
   	     params.dir  = 0;
   	     params.i_s  = 0;
   	     params.k_s  = 0;
   	     params.N    = mesh.x+2;
   	     MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);
  	     Kernel_Launcher(d_CellArray,d_params,dimGridX.x,dimGridX.y,dimGridX.z,
  	                     dimBlock.x,dimBlock.y,dimBlock.z,16000,h_CurrentPeriodic,
   	                     "CurrentPeriodic");

//	     GPU_CurrentPeriodic<<<dimGridX,dimBlock>>>(d_CellArray,d_Jz,2,0,0,0,mesh.x+2);
	     err = getLastError();
	     ThreadSynchronize();
//	     TEST_ERROR(err);
	     checkGPUArray(Jz, d_Jz);
	     PeriodicCurrentBoundaries(Jz,dbgJz,2,1,0,mesh.x+1,0,mesh.z()+1);

	     params.E    = d_Jz;
   	     params.dirE = 2;
   	     params.dir  = 1;
   	     params.i_s  = 0;
   	     params.k_s  = 0;
   	     params.N    = mesh.y+2;
   	     MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);
  	     Kernel_Launcher(d_CellArray,d_params,dimGridY.x,dimGridY.y,dimGridY.z,
  	                     dimBlock.x,dimBlock.y,dimBlock.z,16000,h_CurrentPeriodic,
   	                     "CurrentPeriodic");


//	     GPU_CurrentPeriodic<<<dimGridY,dimBlock>>>(d_CellArray,d_Jz,2,1,0,0,mesh.y+2);
	     err = getLastError();
	     ThreadSynchronize();
//	     TEST_ERROR(err);
	     checkGPUArray(Jz, d_Jz);
	     PeriodicCurrentBoundaries(Jz,dbgJz,2,2,0,mesh.x+1,0,mesh.y+1);

	     params.E    = d_Jz;
   	     params.dirE = 2;
   	     params.dir  = 2;
   	     params.i_s  = 0;
   	     params.k_s  = 0;
   	     params.N    = mesh.y+2;
   	     MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);
  	     Kernel_Launcher(d_CellArray,d_params,dimGridZ.x,dimGridZ.y,dimGridZ.z,
  	                     dimBlock.x,dimBlock.y,dimBlock.z,16000,h_CurrentPeriodic,
   	                     "CurrentPeriodic");


//	     GPU_CurrentPeriodic<<<dimGridZ,dimBlock>>>(d_CellArray,d_Jz,2,2,0,0,mesh.dimz2());
	     err = getLastError();
	     ThreadSynchronize();
//	     TEST_ERROR(err);
	     checkGPUArray(Jz, d_Jz);

	   }

	  void InitQdebug(char *fnjx,char *fnjy,char *fnjz)
	  {

	     read3Darray(fnjx, dbg_Qx);
	     read3Darray(fnjy, dbg_Qy);
	     read3Darray(fnjz, dbg_Qz);
	  }

	  int writeParamsFile(double tex0,double tey0,double tez0,
                          double Tb,double rimp,
                          double rbd,double ni,
	                      double lx,double ly,double lz,
	                      int lp,int nx,int ny,int nz,
	                      double tau,double B0,
	                      double bx,double by,double bz,
	                      double py,double pz,
	                      int beam_plasma,int start_from_file,
	                      int ts,int ms,int phase
	                      )
	  {
		  FILE *f;

		  if((f = fopen("000_params.dat","wt")) == NULL) return 1;

//		  fprintf(f,"\n");
		  fprintf(f,"%15.5e plasma electron temperature along X \n ",tex0);
		  fprintf(f,"%15.5e plasma electron temperature along Y \n ",tey0);
		  fprintf(f,"%15.5e plasma electron temperature along Z \n ",tez0);
		  fprintf(f,"%15.5e beam impulse\n ",rimp);
		  fprintf(f,"%15.5e beam velocity dispersion \n ",Tb);
		  fprintf(f,"%15.5e beam and plasma density ratio \n ",rbd);
		  fprintf(f,"%15.5e plasma density \n ",ni);
		  fprintf(f,"%15.5e external magnetic field (along X) \n ",B0);
		  fprintf(f,"%15.5e domain size X \n ",lx);
		  fprintf(f,"%15.5e domain size Y \n ",ly);
		  fprintf(f,"%15.5e domain size Z \n ",lz);
		  fprintf(f,"%15.5e plasma size Y \n ",py);
		  fprintf(f,"%15.5e plasma size Z \n ",pz);
		  fprintf(f,"%15.5e beam size X \n ",bx);
		  fprintf(f,"%15.5e beam size Y \n ",by);
		  fprintf(f,"%15.5e beam size Z \n ",bz);
		  fprintf(f,"%15d   average number of particles in cell \n ",lp);
		  fprintf(f,"%15d   number of mesh nodes along X \n ",nx);
		  fprintf(f,"%15d   number of mesh nodes along Y \n ",ny);
		  fprintf(f,"%15d   number of mesh nodes along Z \n ",nz);
		  fprintf(f,"%15.5e timestep \n ",tau);
		  fprintf(f,"%15d   1 if beam-plasma interaction, 0 if beam-beam \n",beam_plasma);
		  fprintf(f,"%15d   moment to start from saved\n ",start_from_file);
		  fprintf(f,"%15d   phase to start from save\n ",phase);
		  fprintf(f,"%15d   total steps \n ",ts);
		  fprintf(f,"%15d   number of steps between diagnostic files \n ",ms);


          fclose(f);
          return 0;
	  }


	  int ClearPlasmaBounds(double plasma_ly,double plasma_lz,
			  double *xi,double *yi, double *zi,double *ui,double *vi, double *wi,int *jm)
	  {
		  double plasma_y_max,plasma_y_min,plasma_sh,plasma_z_max,plasma_z_min,plasma_shz,ly,lz;

		   ly = xmax.y;
		   lz = xmax.z();

		   plasma_sh = (ly - plasma_ly)/2;
		   plasma_y_max = ly - plasma_sh;
		   plasma_y_min = plasma_sh;

		   plasma_shz   = (lz - plasma_lz)/2;
		   plasma_z_max = lz - plasma_shz;
		   plasma_z_min = plasma_shz;

		  for(int n = 0;n < *jm;n++)
		  {
			  double y = yi[n];
			  double z = zi[n];

			  if((plasma_y_min > y) || (plasma_y_max < y) || (plasma_z_min > z) || (plasma_z_max < z))
			  {
				  xi[n] = xi[*jm-1];
				  yi[n] = yi[*jm-1];
				  zi[n] = zi[*jm-1];
				  ui[n] = ui[*jm-1];
				  vi[n] = vi[*jm-1];
				  wi[n] = wi[*jm-1];
				  (*jm)--;

				  n--;
			  }
		  }
		  return 0;
	  }

	  void LoadTestData(int nt,int part_nt)
	  {
	     vector<Particle> vp,bin_vp;
	     char d_exfile[100],d_eyfile[100],d_ezfile[100],d_hxfile[100],d_hyfile[100],d_hzfile[100];
	     char d_0exfile[100],d_0eyfile[100],d_0ezfile[100];
	     char jxfile[100],jyfile[100],jzfile[100];
	     char np_jxfile[100],np_jyfile[100],np_jzfile[100];
	     char np_exfile[100],np_eyfile[100],np_ezfile[100];
	     char d_jxfile[100],d_jyfile[100],d_jzfile[100];
	     char qxfile[100],qyfile[100],qzfile[100];
	     char pfile[100],nextpfile[100];
	     char part_name[100];

	     sprintf(qxfile,"dnqx%06d.dat",nt);
	     sprintf(qyfile,"dnqy%06d.dat",nt);
	     sprintf(qzfile,"dnqz%06d.dat",nt);

	     readDebugArray("hxlg",Hx,nt,5);
	     readDebugArray("hylg",Hy,nt,5);
	     readDebugArray("hzlg",Hz,nt,5);

	     sprintf(d_exfile,"dnex%06d.dat",2*nt-1);
	     sprintf(d_eyfile,"dney%06d.dat",2*nt-1);
	     sprintf(d_ezfile,"dnez%06d.dat",2*nt-1);

	     sprintf(d_0exfile,"dnex%06d.dat",2*nt-2);
	     sprintf(d_0eyfile,"dney%06d.dat",2*nt-2);
	     sprintf(d_0ezfile,"dnez%06d.dat",2*nt-2);

	     sprintf(d_hxfile,"dnhx%06d.dat",2*nt-1);
	     sprintf(d_hyfile,"dnhy%06d.dat",2*nt-1);
	     printf("%s \n",d_hyfile);
	     sprintf(d_hzfile,"dnhz%06d.dat",2*nt-1);

	     sprintf(jxfile,"dnjx%06d.dat",2*nt);
	     sprintf(jyfile,"dnjy%06d.dat",2*nt);
	     sprintf(jzfile,"dnjz%06d.dat",2*nt);

	     sprintf(d_jxfile,"npjx%06d.dat",2*nt);
	     sprintf(d_jyfile,"npjy%06d.dat",2*nt);
	     sprintf(d_jzfile,"npjz%06d.dat",2*nt);

	     sprintf(np_jxfile,"npjx%06d.dat",2*nt);
	     sprintf(np_jyfile,"npjy%06d.dat",2*nt);
	     sprintf(np_jzfile,"npjz%06d.dat",2*nt);

	     sprintf(np_exfile,"exlg%03d.dat",2*nt);
	     sprintf(np_eyfile,"eylg%03d.dat",2*nt);
	     sprintf(np_ezfile,"ezlg%03d.dat",2*nt);

	     sprintf(pfile,    "part%06d000.dat",nt);
	     sprintf(nextpfile,"part%06d000.dat",nt+2);

	     InitQdebug(qxfile,qyfile,qzfile);

	     InitCurrents(jxfile,jyfile,jzfile,d_jxfile,d_jyfile,d_jzfile,np_jxfile,np_jyfile,np_jzfile,0);

	     if(nt > 1)
	     {
	    	 ClearAllParticles();

	     }

	     sprintf(part_name,"mumu%03d%05d.dat",start_phase,nt_start_from_file);

   	     double *dbg_ion_x, *dbg_ion_y, *dbg_ion_z,
   	            *dbg_ion_px,*dbg_ion_py,*dbg_ion_pz,ion_q_m,ion_m;
   	     int total_ions;
   	     double *dbg_el_x, *dbg_el_y, *dbg_el_z,
   	            *dbg_el_px,*dbg_el_py,*dbg_el_pz,el_q_m,el_m;
   	     int total_electrons;
   	     double *dbg_beam_x, *dbg_beam_y, *dbg_beam_z,
   	            *dbg_beam_px,*dbg_beam_py,*dbg_beam_pz,beam_q_m,beam_m;
	     int total_beam_electrons;
	     //int total = 160000;

//	     tex0 = 1e-3;
//	     tey0 = 1e-3;
//	     tez0 = 1e-3;
	     tol  = 1e-15;

//	     Tb   = 0.0;
//	     rimp = 0.2;
//	     rbd  = 1.0e-2;
//	     ni   = 1.0;
//	     meh  = 0;
//	     lp   = 100;



	     ParticleArrays ions,electrons,beam_electrons;

/////////////////////////////////////////////////////////////////////
	  if(nt_start_from_file >= 0)
	  {
	     InitBinaryParticlesArrays(part_name,part_nt,
	    		 &ions,&electrons,&beam_electrons,
	    		 mesh.x,mesh.y,mesh.z(),
	    		 beam_plasma);
	     change_work_dir();
	  }

	  if(beam_plasma == 1)
	  {

	     initial[0].total    = N;
	     initial[1].total    = 2*N;
	     initial[2].total    = N;

	     diagnostics[0].total    = N;
	     diagnostics[1].total    = 2*N;
	     diagnostics[2].total    = N;

	     sorts = 3;
	  }
	  else
	  {
		  if(nt_start_from_file >= 0)
		  {
		     initial[0].total = ions.total;
		  }
		  else
		  {
			 initial[0].total = N;
		  }
		  sorts = 1;
	  }




	     component_total = 2*N;

//	     diagnose = (float *)malloc(sizeof(float)*sorts*component_total*8);
//	     MemoryAllocate((void**)&d_diagnose,sizeof(float)*sorts*component_total*8);

	     AllocateBinaryParticlesArrays(&(initial[0]),&(initial[1]),&(initial[2]));
	     AllocateBinaryParticlesArraysFloat(&(diagnostics[0]),&(diagnostics[1]),&(diagnostics[2]));

//	     AllocateDeviceParticleDiagnosticPointers(d_diagnostics,diagnostics);

//	  for(int j = 535;j < 542;j++)
//	  {
//	     printf("%10d %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e \n",j,
//	    		 electrons.dbg_x[j],electrons.dbg_y[j],electrons.dbg_z[j],
//   	    		 electrons.dbg_px[j],electrons.dbg_py[j],electrons.dbg_pz[j]
//	   		                                           );
//	  }
         writeParamsFile(tex0,tey0,tez0,
        		         Tb,rimp,rbd,ni,
                         xmax.x,xmax.y,xmax.z(),
	                     lp,mesh.x,mesh.y,mesh.z(),
	                     tau,Bx0,beam_max.x,beam_max.y,beam_max.z(),
	                     plasma_dim_y,plasma_dim_z,
	                     beam_plasma,nt_start_from_file,
	                     total_steps,minor_steps,start_phase);



      int jmb_real;

      if((nt_start_from_file < 0) && (beam_plasma == 1))
      {
	     InitUniformMaxwellianParticles(1,N,tex0,tey0,tez0,
	    		   beam_max.x,beam_max.y,beam_max.z(),&jmb_real,
	    		   xmax.x,xmax.y,xmax.z(),meh,Tb,rimp,rbd,
				   initial[0].dbg_x,initial[0].dbg_y,initial[0].dbg_z,
				   initial[0].dbg_px,initial[0].dbg_py,initial[0].dbg_pz,
 				   initial[2].dbg_x,initial[2].dbg_y,initial[2].dbg_z,
 				   initial[2].dbg_px,initial[2].dbg_py,initial[2].dbg_pz,
 				   initial[1].dbg_x,initial[1].dbg_y,initial[1].dbg_z,
 				   initial[1].dbg_px,initial[1].dbg_py,initial[1].dbg_pz);

	     ClearPlasmaBounds(plasma_dim_y,plasma_dim_z,
	    		   initial[0].dbg_x, initial[0].dbg_y, initial[0].dbg_z,
				   initial[0].dbg_px,initial[0].dbg_py,initial[0].dbg_pz,&(initial[0].total));
	     ClearPlasmaBounds(plasma_dim_y,plasma_dim_z,
	    		   initial[1].dbg_x, initial[1].dbg_y, initial[1].dbg_z,
				   initial[1].dbg_px,initial[1].dbg_py,initial[1].dbg_pz,&(initial[1].total));
	     ClearPlasmaBounds(plasma_dim_y,plasma_dim_z,
	    		   initial[2].dbg_x, initial[2].dbg_y, initial[2].dbg_z,
				   initial[2].dbg_px,initial[2].dbg_py,initial[2].dbg_pz,&jmb_real);
      }


#ifdef DEBUG_LOAD_PARTICLES_FROM_FILE
	       double t_beam_electron_x = compare(beam_electrons1.dbg_x,beam_electrons.dbg_x,N,"beam_electronx",tol);
	       double t_beam_electron_y = compare(beam_electrons1.dbg_y,beam_electrons.dbg_y,N,"beam_electrony",tol);
	       double t_beam_electron_z = compare(beam_electrons1.dbg_z,beam_electrons.dbg_z,N,"beam_electronz",tol);
	       double t_beam_electron_px = compare_prints(beam_electrons1.dbg_px,beam_electrons.dbg_px,10,"beam_electronpx",tol,1);
	       double t_beam_electron_py = compare(beam_electrons1.dbg_py,beam_electrons.dbg_py,N,"beam_electronpy",tol);
	       double t_beam_electron_pz = compare(beam_electrons1.dbg_pz,beam_electrons.dbg_pz,N,"beam_electronpz",tol);
	       printf("ELECTRONS x %15.6e y %15.6e z %15.6e px %15.6e py %15.6e pz %15.6e \n",t_beam_electron_x,t_beam_electron_y,t_beam_electron_z,t_beam_electron_px,
	     	 t_beam_electron_py,t_beam_electron_pz);
#endif
//	       double t_electron_x = compare(electrons1.dbg_x,electrons.dbg_x,N,"electronx",tol);
//	       double t_electron_y = compare(electrons1.dbg_y,electrons.dbg_y,N,"electrony",tol);
//	       double t_electron_z = compare(electrons1.dbg_z,electrons.dbg_z,N,"electronz",tol);
//	       double t_electron_px = compare_prints(electrons1.dbg_px,electrons.dbg_px,10,"electronpx",tol,1);
//	       double t_electron_py = compare(electrons1.dbg_py,electrons.dbg_py,N,"electronpy",tol);
//	       double t_electron_pz = compare(electrons1.dbg_pz,electrons.dbg_pz,N,"electronpz",tol);
//	       printf("BEAM x %15.6e y %15.6e z %15.6e px %15.6e py %15.6e pz %15.6e \n",t_electron_x,t_electron_y,t_electron_z,t_electron_px,
//	     	 t_electron_py,t_electron_pz);
         if(beam_plasma == 1)
         {
           getMassCharge(&initial[0],&initial[1],&initial[2],ni,rbd,lp);


  	       diagnostics[0].m        = initial[0].m[0];
  	       diagnostics[1].m        = initial[1].m[0];
  	       diagnostics[2].m        = initial[2].m[0];
  	       diagnostics[0].q_m      = initial[0].q_m;
  	       diagnostics[1].q_m      = initial[1].q_m;
  	       diagnostics[2].q_m      = initial[2].q_m;

	       AllocateDeviceParticleDiagnosticPointers(&d_diagnostics,
	     			                              &host_copy_d_diagnostics,
	     			                              &diagnostics);

	       initial[2].total = jmb_real;
	       diagnostics[2].total = jmb_real;
         }
///////////////////////////////////////////////////////////////////////
         if(beam_plasma == 1)
         {
	        InitBinaryParticles(bin_vp,part_nt,
	     			  initial[0].dbg_x, initial[0].dbg_y, initial[0].dbg_z,
	     			  initial[0].dbg_px,initial[0].dbg_py,initial[0].dbg_pz,
	     			  initial[0].q_m,initial[0].m[0],
	     			  initial[0].total,
	     			  initial[1].dbg_x,initial[1].dbg_y,initial[1].dbg_z,
	     			  initial[1].dbg_px,initial[1].dbg_py,initial[1].dbg_pz,
	     			  initial[1].q_m,initial[1].m[0],
	     			  initial[1].total,
	     			  initial[2].dbg_x,initial[2].dbg_y,initial[2].dbg_z,
	     			  initial[2].dbg_px,initial[2].dbg_py,initial[2].dbg_pz,
	     			  initial[2].q_m,initial[2].m[0],
	     			  initial[2].total);
         }
         else
         {
 	        InitBinaryParticlesTotal(bin_vp,part_nt,
 	     			  initial[0].dbg_x, initial[0].dbg_y, initial[0].dbg_z,
 	     			  initial[0].dbg_px,initial[0].dbg_py,initial[0].dbg_pz,
 	     			  initial[0].q_m,initial[0].m,
 	     			  initial[0].total);
         }
//	      InitBinaryParticles(part_name,bin_vp,part_nt);


         AssignArraysToCells();

	  }

void readParticles(char *pfile,char *nextpfile)
{
	vector<Particle> vp;

	if(!strncmp(pfile,"mumu",4))
	{
//		InitBinaryParticles(pfile,vp,4);
	}
	else
	{
       InitParticles(pfile,vp);
	}
#ifdef DEBUG_PLASMA
	if(!strncmp(pfile,"mumu",4))
	{
	    InitBinaryParticlesNext(nextpfile, vp);
	}
	else
	{
	    InitParticlesNext(nextpfile, vp);
	}
#endif
    Distribute(vp);

    AssignArraysToCells();
}

void AssignCellsToArraysGPU()
{
	int err;
	dim3 dimGrid(mesh.x,mesh.y,mesh.z()),dimBlockExt(CellExtent,CellExtent,CellExtent);




	params.d_Ex = d_Ex;
    params.d_Ey = d_Ey;
	params.d_Ez = d_Ez;

	params.d_Hx = d_Hx;
    params.d_Hy = d_Hy;
	params.d_Hz = d_Hz;
	//		static int first  = 1;
    MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

	Kernel_Launcher(d_CellArray,d_params,mesh.x,mesh.y,mesh.z(),
			CellExtent,CellExtent,CellExtent,16000,h_GPU_SetFieldsToCells_SingleNode,
	               "SetFieldsToCells");




//	GPU_SetFieldsToCells<<<dimGrid, dimBlockExt>>>(d_CellArray,d_Ex,d_Ey,d_Ez,d_Hx,d_Hy,d_Hz);
    err = getLastError();
    ThreadSynchronize();
//    TEST_ERROR(err);

}


	  void AssignCellsToArrays()
	{

	     for(int n = 0;n < (*AllCells).size();n++)
	     {
	         Cell<Particle,dims>  c = (*AllCells)[n];
		     c.writeAllToArrays(Jx,Jy,Jz,Rho,0);
	     }
	     CheckArray(Jx, dbgJx);
	     SetPeriodicCurrents(0);
	     CheckArray(Jx, dbgJx);
	}
	  void AssignArraysToCells()
	  {
	     for(int n = 0;n < (*AllCells).size();n++)
	     {

	         Cell<Particle,dims> c = (*AllCells)[n];
		     c.readFieldsFromArrays(Ex,Ey,Ez,Hx,Hy,Hz);
	     }
	  }

	  void ParticleLog()
	{
	#ifndef DEBUG_PLASMA
	     return;
	#endif

	     FILE *f;
	     char  fname[100];
	     int   num = 0;

	     sprintf(fname,"particles.dat");

	     if((f = fopen(fname,"wt")) == NULL) return;

	     for(int n = 0;n < (*AllCells).size();n++)
	     {
	         Cell<Particle,dims>  c = (*AllCells)[n];
	#ifdef GPU_PARTICLE
	   	 vector<Particle>  pvec_device;// = c.GetParticles();
	   	 vector<Particle> pvec = pvec_device;
	#else
		 vector<Particle>  pvec = c.GetParticles();
	#endif

		 for(int i = 0;i < pvec.size();i++)
		 {
		      Particle p = pvec[i];

		      p.Print(f,num++);
		 }

	     }

	     fclose(f);
	  }


	  void write3Darray(char *name,double *d)
	  {
	    char fname[100];
	    GPUCell<Particle,dims> c = (*AllCells)[0];
	    FILE *f;

	    sprintf(fname,"%s_fiel3d.dat",name);

	    if((f = fopen(fname,"wt")) == NULL) return;

	    for(int i = 1;i < mesh.x+1;i++)
	    {
	        for(int l = 1;l < mesh.y+1;l++)
	        {
	            for(int k = 1;k < mesh.z()+1;k++)
		    {
		        int n = c.getGlobalCellNumber(i,l,k);

			fprintf(f,"%15.5e %15.5e %15.5e %25.15e \n",c.getNodeX(i),c.getNodeY(l),c.getNodeZ(k),d[n]);
		    }
		}
	    }

	    fclose(f);
	}

void write3D_GPUArray(char *name,double *d_d)
{
	double *d;

#ifndef WRITE_3D_DEBUG_ARRAYS
	return;
#endif

	d = (double *)malloc(sizeof(double)*mesh.size2());

	int err = MemoryCopy(d,d_d,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);

	write3Darray(name,d);
}

void readControlPoint(FILE **f1,char *fncpy,int num,int nt,int part_read,int field_assign,
		double *ex,double *ey,double *ez,
		double *hx,double *hy,double *hz,
		double *jx,double *jy,double *jz,
		double *qx,double *qy,double *qz,
		double *x,double *y,double *z,
		double *px,double *py,double *pz
		)
{
	char fn[100],fn_next[100];
	FILE *f;

	sprintf(fn,"mumu%03d%08d.dat",num,nt);
	strcpy(fncpy,fn);
	sprintf(fn_next,"mumu%03d%05d.dat",num,nt+1);
	if((f = fopen(fn,"rb")) == NULL) return;
	if(part_read)
	{
	   *f1 = f;
	}

	readFortranBinaryArray(f,ex);
	readFortranBinaryArray(f,ey);
	readFortranBinaryArray(f,ez);
	readFortranBinaryArray(f,hx);
	readFortranBinaryArray(f,hy);
	readFortranBinaryArray(f,hz);
	readFortranBinaryArray(f,jx);
	readFortranBinaryArray(f,jy);
	readFortranBinaryArray(f,jz);

	readFortranBinaryArray(f,qx);
	readFortranBinaryArray(f,qy);
	readFortranBinaryArray(f,qz);

	if(field_assign == 1) AssignArraysToCells();

}

double checkControlMatrix(string wh,int nt,string name, double *d_m)
{
	double /*t_ex,t_ey,t_ez,t_hx,t_hy,t_hz,*/t_jx,t_jy,t_jz;
	char fn[100];//,fn_next[100];
	FILE *f;

#ifndef CHECK_CONTROL_MATRIX
	return 0.0;
#endif

	sprintf(fn,"wcmx_%4s_%08d_%2s.dat",wh.c_str(),nt,name.c_str());
	if((f = fopen(fn,"rb")) == NULL) return -1.0;

	readFortranBinaryArray(f,dbgJx);

	t_jx = checkGPUArray(dbgJx,d_m);

    return t_jx;
}


void checkCurrentControlPoint(int j,int nt)
{
	 double /*t_ex,t_ey,t_ez,t_hx,t_hy,t_hz,*/t_jx,t_jy,t_jz;
		char fn[100];//,fn_next[100];
		FILE *f;

		sprintf(fn,"curr%03d%05d.dat",nt,j);
		if((f = fopen(fn,"rb")) == NULL) return;

		readFortranBinaryArray(f,dbgJx);
		readFortranBinaryArray(f,dbgJy);
		readFortranBinaryArray(f,dbgJz);

	 t_jx = CheckArraySilent(Jx,dbgJx);
	 t_jy = CheckArraySilent(Jy,dbgJy);
	 t_jz = CheckArraySilent(Jz,dbgJz);

     printf("Jx %15.5e Jy %15.5e Jz %15.5e \n",t_jx,t_jy,t_jz);
}

void checkControlPoint(int num,int nt,int check_part)
{
	 double t_ex,t_ey,t_ez,t_hx,t_hy,t_hz,t_jx,t_jy,t_jz;
	 double t_qx,t_qy,t_qz,t_njx,t_njy,t_njz;
	 FILE *f;
	 char fn_copy[100];
	 struct sysinfo info;

	 if(nt == 0 && num > 0) return;

	 double cf = checkFields(d_Ex,d_Ey,d_Ez,d_Hx,d_Hy,d_Hz,
		 			d_Jx,d_Jy,d_Jz,d_Qx,d_Qy,d_Qz,
		 		    1,mesh.size2(),num,nt);

	 memory_monitor("checkControlPoint1",nt);


	 if(nt % FORTRAN_NUMBER_OF_SMALL_STEPS != 0) return;

	 memory_monitor("checkControlPoint2",nt);

	 readControlPoint(&f,fn_copy,num,nt,1,0,dbgEx,dbgEy,dbgEz,dbgHx,dbgHy,dbgHz,dbgJx,dbgJy,dbgJz,dbg_Qx,dbg_Qy,dbg_Qz,
                     dbg_x,dbg_y,dbg_z,dbg_px,dbg_py,dbg_pz);

	 memory_monitor("checkControlPoint3",nt);


	 t_ex = CheckArraySilent(Ex,dbgEx);
	 t_ey = CheckArraySilent(Ey,dbgEy);
	 t_ez = CheckArraySilent(Ez,dbgEz);
	 t_hx = CheckArraySilent(Hx,dbgHx);
	 t_hy = CheckArraySilent(Hy,dbgHy);
	 t_hz = CheckArraySilent(Hz,dbgHz);
	 t_jx = CheckArraySilent(Jx,dbgJx);
	 t_jy = CheckArraySilent(Jy,dbgJy);
	 t_jz = CheckArraySilent(Jz,dbgJz);

	 memory_monitor("checkControlPoint4",nt);

	 t_ex = CheckGPUArraySilent(dbgEx,d_Ex);
	 t_ey = CheckGPUArraySilent(dbgEy,d_Ey);
	 t_ez = CheckGPUArraySilent(dbgEz,d_Ez);
	 t_hx = CheckGPUArraySilent(dbgHx,d_Hx);
	 t_hy = CheckGPUArraySilent(dbgHy,d_Hy);
	 t_hz = CheckGPUArraySilent(dbgHz,d_Hz);

	 t_qx = CheckGPUArraySilent(dbg_Qx,d_Qx);
	 t_qy = CheckGPUArraySilent(dbg_Qy,d_Qy);
	 t_qz = CheckGPUArraySilent(dbg_Qz,d_Qz);

	 t_jx = CheckGPUArraySilent(dbgJx,d_Jx);
	 t_jy = CheckGPUArraySilent(dbgJy,d_Jy);
	 t_jz = CheckGPUArraySilent(dbgJz,d_Jz);

	 t_njx = CheckGPUArraySilent(dbgJx,d_Jx);
	 t_njy = CheckGPUArraySilent(dbgJy,d_Jy);
	 t_njz = CheckGPUArraySilent(dbgJz,d_Jz);

	 memory_monitor("checkControlPoint5",nt);

	 double t_cmp_jx = checkGPUArray(dbgJx,d_Jx,"Jx","step",nt);
	 double t_cmp_jy = checkGPUArray(dbgJy,d_Jy,"Jy","step",nt);
	 double t_cmp_jz = checkGPUArray(dbgJz,d_Jz,"Jz","step",nt);

#ifdef CONTROL_DIFF_GPU_PRINTS
     printf("GPU: Ex %15.5e Ey %15.5e Ez %15.5e \n",t_ex,t_ey,t_ez);
     printf("GPU: Hx %15.5e Hy %15.5e Ez %15.5e \n",t_hx,t_hy,t_hz);
     printf("GPU: Jx %15.5e Jy %15.5e Jz %15.5e \n",t_jx,t_jy,t_jz);
     printf("GPU compare : Jx %15.5e Jy %15.5e Jz %15.5e \n",t_cmp_jx,t_cmp_jy,t_cmp_jz);
#endif

     memory_monitor("checkControlPoint6",nt);

     double cp = checkControlPointParticles(num,f,fn_copy,nt);

     f_prec_report = fopen("control_points.dat","at");
     fprintf(f_prec_report,"nt %5d num %3d Ex %15.5e Ey %15.5e Ez %15.5e Hx %15.5e Hy %15.5e Hz %15.5e Jx %15.5e Jy %15.5e Jz %15.5e Qx %15.5e Qy %15.5e Qz %15.5e particles %15.5e\n",
    		 nt,num,
    		 t_ex,t_ey,t_ez,
    		 t_hx,t_hy,t_hz,
    		 t_jx,t_jy,t_jz,
    		 t_qx,t_qy,t_qz,
    		 cp
    		 );
     fclose(f_prec_report);

     memory_monitor("checkControlPoint7",nt);

     fclose(f);
}

//	int readFortranBinaryArray(FILE *f, double* d)
//	{
////	    Cell<Particle,dims>  c = (*AllCells)[0];
//	    int t,err;//,n;
//	    fread(&t,sizeof(int),1,f);
//	     if((err = ferror(f)) != 0)
//	    	 {
//	    	 	 return err ;
//	    	 }
//
//	    fread(d,1,t,f);
//	     if((err = ferror(f)) != 0)
//	    	 {
//	    	 	 return err ;
//	    	 }
//
//	    fread(&t,sizeof(int),1,f);
//	     if((err = ferror(f)) != 0)
//	    	 {
//	    	 	 return err ;
//	    	 }
//
//
//
//#ifdef READ_DEBUG_PRINTS
//	    for(int i = 1; i <= mesh.x+2;i++)
//	    {
//	    	for(int l = 1; l <= mesh.y+2;l++)
//	    	{
//	    		for(int k = 1;k <= mesh.dimz2();k++)
//	    		{
//	    			n = c.getFortranCellNumber(i,l,k);
//	    			printf("%5d %5d %5d %25.15e \n",i,l,k,d[n]);
//	    		}
//	    	}
//	    }
//#endif
//
//	    	    return t;
//	}


	  void read3Darray(char* name, double* d)
	{
	    char str[100];
	    Cell<Particle,dims>  c = (*AllCells)[0];
	    FILE *f;

	    if((f = fopen(name,"rt")) == NULL) return;

	    while(fgets(str,100,f) != NULL)
	    {
	          int i = atoi(str);
	          int l = atoi(str + 10);
	          int k = atoi(str + 20);
		  double t = atof(str + 30);
		  int i1,l1,k1,n = c.getFortranCellNumber(i,l,k);
		  c.getFortranCellTriplet(n,&i1,&l1,&k1);
		  d[n] = t;
	    }

	    fclose(f);

	}

	int readDebugArray(string name, double* d,int nt,int col)
	{
		char dfile[100];

		if(!strncmp(name.c_str()+2,"lg",2))
		{
			sprintf(dfile,"%s%03d.dat",name.c_str(),nt);
			if(name[0] == 'e')
			{
			   read3DarrayLog(dfile, d,50,col);
			}
			else
			{
				read3DarrayLog(dfile,d,40,col);
			}
		}
		else
		{
			sprintf(dfile,"%s%06d.dat",name.c_str(),nt);
			read3Darray(dfile,d);
		}
		return 0;
	}


	  void read3DarrayModified(char* name, double* d,double a)
	  {
	      char str[100];
	      Cell<Particle,dims>  c = (*AllCells)[0];
	      FILE *f;

	      //sprintf(fname,"%s_fiel3d.dat",name);

	      if((f = fopen(name,"rt")) == NULL) return;

	      while(fgets(str,100,f) != NULL)
	      {
	            int i = atoi(str) - 1;
	            int l = atoi(str + 10) - 1;
	            int k = atoi(str + 20) - 1;
	  	  double t = atof(str + 30);
	  	  int n = c.getGlobalCellNumber(i,l,k);
	  	   t = a*(i + 100*l + 10000*k);
	  	   d[n] = t;
	      }

	      fclose(f);

	  }



	  void write3DcellArray(char *name,int code)
	  {
	    char fname[100];
	    Cell<Particle,dims> & c = AllCells[0];
	    FILE *f;
#ifndef WRITE_3D_CELL_ARRAY
	    return;
#endif

	    sprintf(fname,"%03d_cells.dat",code,name);

	    if((f = fopen(fname,"wt")) == NULL) return;

	    for(int i = 1;i < mesh.x+1;i++)
	    {
	        for(int l = 1;l < mesh.y+1;l++)
	        {
	            for(int k = 1;k < mesh.z()+1;k++)
		    {
		        int n = c.getGlobalCellNumber(i,l,k);
			Cell<Particle,dims> & cc = AllCells[n];

			fprintf(f,"%15.5e %15.5e %15.5e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",c.getNodeX(i),c.getNodeY(l),c.getNodeZ(k),
				cc.getCoreCell(code,1,1,1),
				cc.getCoreCell(code,1,1,2),
				cc.getCoreCell(code,1,2,1),
				cc.getCoreCell(code,1,2,2),
				cc.getCoreCell(code,2,1,1),
				cc.getCoreCell(code,2,1,2),
				cc.getCoreCell(code,2,2,1),
				cc.getCoreCell(code,2,2,2)
			);
		    }
		}
	    }

	    fclose(f);
	}



void copyCellCurrentsToDevice(CellDouble *d_jx,CellDouble *d_jy,CellDouble *d_jz,
		                      CellDouble *h_jx,CellDouble *h_jy,CellDouble *h_jz)
{
	int err;

 	err = MemoryCopy(d_jx,h_jx,sizeof(CellDouble),HOST_TO_DEVICE);
 	if(err != 0)
 	        {
 	         	printf("1copyCellCurrentsToDevice err %d %s \n",err,getErrorString(err));
 	       	exit(0);
 	        }
 	err = MemoryCopy(d_jy,h_jy,sizeof(CellDouble),HOST_TO_DEVICE);
 	if(err != 0)
 	        {
 	         	printf("2copyCellCurrentsToDevice err %d %s \n",err,getErrorString(err));
 	       	exit(0);
 	        }
 	err = MemoryCopy(d_jz,h_jz,sizeof(CellDouble),HOST_TO_DEVICE);
 	if(err != 0)
 	        {
 	         	printf("3copyCellCurrentsToDevice err %d %s \n",err,getErrorString(err));
 	       	exit(0);
 	        }

}


#include "check_field.h"


//double CheckArray	(double* a, double* dbg_a,FILE *f)
//	{
//	    Cell<Particle,dims> c = (*AllCells)[0];
//	    int wrong = 0;
//	    double diff = 0.0;
//
//
//
////#ifdef CHECK_ARRAY_DETAIL_PRINTS
//	    fprintf(f,"begin array checking=============================\n");
////#endif
//	    for(int n = 0;n < mesh.size2();n++)
//	    {
////	        double t  = a[n];
////		    double dt = dbg_a[n];
//            diff += pow(a[n] - dbg_a[n],2.0);
//
//	        if(fabs(a[n] - dbg_a[n]) > TOLERANCE)
//		    {
//
//		       int3 i = c.getCellTripletNumber(n);
//#ifdef CHECK_ARRAY_DETAIL_PRINTS
//		       fprintf(f,"n %5d i %3d l %3d k %3d %15.5e dbg %15.5e diff %15.5e wrong %10d \n",
//				   n,i.x+1,i.y+1,i.z+1,a[n],dbg_a[n],fabs(a[n] - dbg_a[n]),wrong++);
//#endif
//     		}
//	    }
//#ifdef CHECK_ARRAY_DETAIL_PRINTS
//	    fprintf(f,"  end array checking============================= %.4f less than %15.5e diff %15.5e \n",
//	    		(1.0-((double)wrong/(mesh.size2()))),TOLERANCE,
//	    		pow(diff/(mesh.size2()),0.5)
//	    	  );
//#endif
//
//	    return (1.0-((double)wrong/(mesh.size2())));
//	}
//
//double CheckArray	(double* a, double* dbg_a)
//	{
//	    Cell<Particle,dims> c = (*AllCells)[0];
//	    int wrong = 0;
//	    double diff = 0.0;
//#ifdef CHECK_ARRAY_DETAIL_PRINTS
//	    puts("begin array checking2=============================");
//#endif
//	    for(int n = 0;n < mesh.size2();n++)
//	    {
////	        double t  = a[n];
////		    double dt = dbg_a[n];
//            diff += pow(a[n] - dbg_a[n],2.0);
//
//	        if(fabs(a[n] - dbg_a[n]) > TOLERANCE)
//		    {
//
//		       int3 i = c.getCellTripletNumber(n);
//#ifdef CHECK_ARRAY_DETAIL_PRINTS
//		       printf("n %5d i %3d l %3d k %3d %15.5e dbg %15.5e diff %15.5e wrong %10d \n",
//				   n,i.x+1,i.y+1,i.z+1,a[n],dbg_a[n],fabs(a[n] - dbg_a[n]),wrong++);
//#endif
//     		}
//	    }
//#ifdef CHECK_ARRAY_DETAIL_PRINTS
//	    printf("  end array checking============================= %.4f less than %15.5e diff %15.5e \n",
//	    		(1.0-((double)wrong/(mesh.size2()))),TOLERANCE,
//	    		pow(diff/(mesh.size2()),0.5)
//	    	  );
//#endif
//
//	    return (1.0-((double)wrong/(mesh.size2())));
//	}
//
//
//double CheckArraySilent	(double* a, double* dbg_a)
//	{
//	    Cell<Particle,dims> c = (*AllCells)[0];
////	    int wrong = 0;
//	    double diff = 0.0;
//	   // puts("begin array checking=============================");
//	    for(int n = 0;n < mesh.size2();n++)
//	    {
////	        double t  = a[n];
////		    double dt = dbg_a[n];
//            diff += pow(a[n] - dbg_a[n],2.0);
//
//	        if(fabs(a[n] - dbg_a[n]) > TOLERANCE)
//		    {
//
//		       int3 i = c.getCellTripletNumber(n);
//
////		       printf("n %5d i %3d l %3d k %3d %15.5e dbg %15.5e diff %15.5e wrong %10d \n",
////				   n,i.x+1,i.y+1,i.z+1,a[n],dbg_a[n],fabs(a[n] - dbg_a[n]),wrong++);
//     		}
//	    }
////	    printf("  end array checking============================= %.2f less than %15.5e diff %15.5e \n",
////	    		(1.0-((double)wrong/((mesh.x + 2)*(mesh.y + 2)*(mesh.z() + 2)))),TOLERANCE,
////	    		pow(diff/((mesh.x + 2)*(mesh.y + 2)*(mesh.z() + 2)),0.5)
////	    	  );
//
//	    return pow(diff/((mesh.x + 2)*(mesh.y + 2)*(mesh.z() + 2)),0.5);
//	}
//
//double CheckGPUArraySilent	(double* a, double* d_a)
//	{
//	    static double *t;
//	    static int f = 1;
//	    cudaError_t err;
//
//	    if(f == 1)
//	    {
//	    	 t = (double *)malloc(sizeof(double)*mesh.size2());
//	    	 f = 0;
//	    }
//	    cudaMemcpy(t,d_a,sizeof(double)*mesh.size2(),cudaMemcpyDeviceToHost);
//	    err = cudaGetLastError();
//	    if(err != cudaSuccess)
//	            {
//	             	printf("CheckArraySilent err %d %s \n",err,cudaGetErrorString(err));
//	            	exit(0);
//	            }
//
//
//	   return CheckArraySilent(a,t);
//	}




	int CheckValue(double *a, double *dbg_a, int n)
	{
	    Cell<Particle,dims>  c = (*AllCells)[0];

	    if(fabs(a[n] - dbg_a[n]) > TOLERANCE)
	    {

	       int3 i = c.getCellTripletNumber(n);
#ifdef CHECK_VALUE_DEBUG_PRINTS
	       printf("value n %5d i %3d l %3d k %3d %15.5e dbg %1.5e diff %15.5e \n",n,i.x,i.y,i.z,a[n],dbg_a[n],fabs(a[n] - dbg_a[n]));
#endif

	       return 0;

	    }

	    return 1;
	}


	void read3DarrayLog(char* name, double* d, int offset, int col)
	{
	    char str[1000];
	    Cell<Particle,dims> c = (*AllCells)[0];
	    FILE *f;

	    //sprintf(fname,"%s_fiel3d.dat",name);

	    if((f = fopen(name,"rt")) == NULL) return;

	    while(fgets(str,1000,f) != NULL)
	    {
	          //str += offset;

	          int i = atoi(str + offset)      - 1;
	          int l = atoi(str + offset + 5)  - 1;
	          int k = atoi(str + offset + 10) - 1;
		  double t = atof(str + offset + 15 + col*25);
		  int n = c.getGlobalCellNumber(i,l,k);
		  d[n] = t;
#ifdef READ_ARRAY_LOG_PRINTS
		  printf("%d %d %5d %5d %15.5e \n",i,l,k,n,t);
#endif
	    }

	    fclose(f);

	}

	void InitParticlesNext(char* fname, vector< Particle >& vp)
	{
	    FILE *f;
	     char str[1000];
	     //double3 x;
	     int n = 0;

	     if((f = fopen(fname,"rt")) == NULL) return;

	     while(fgets(str,1000,f) != NULL)
	     {

	          //px  = atof(str + 75);
	          //py  = atof(str + 100);
	          //pz  = atof(str + 125);
	          //m   = atof(str + 150);
	          //q_m = atof(str + 175);

		  Particle & p = vp[n++];
		  //p.SetXnext(x);
	     }
	}

	void InitBinaryParticlesNext(char *fn, vector< Particle >& vp)
	{
	     FILE *f;
	     double *buf,q_m,m,tp;
	     int n = 0,t;

	     buf = (double *)malloc(sizeof(double)*mesh.size2());

	     if((f = fopen(fn,"rb")) == NULL) return;
	     readFortranBinaryArray(f,buf);
	     readFortranBinaryArray(f,buf);
	     readFortranBinaryArray(f,buf);

	     readFortranBinaryArray(f,buf);
	     readFortranBinaryArray(f,buf);
	     readFortranBinaryArray(f,buf);

	     readFortranBinaryArray(f,buf);
	     	     readFortranBinaryArray(f,buf);
	     	     readFortranBinaryArray(f,buf);

	     	     fread(&t,sizeof(int),1,f);
	     	     fread(&tp,sizeof(double),1,f);
	     	     total_particles = (int)tp;
	     	     fread(&q_m,sizeof(double),1,f);
	     	     fread(&m,sizeof(double),1,f);
	     	     fread(&t,sizeof(int),1,f);
	              dbg_x = (double *)malloc(sizeof(double)*total_particles);
	              dbg_y = (double *)malloc(sizeof(double)*total_particles);
	              dbg_z = (double *)malloc(sizeof(double)*total_particles);
	              dbg_px = (double *)malloc(sizeof(double)*total_particles);
	              dbg_py = (double *)malloc(sizeof(double)*total_particles);
	              dbg_pz = (double *)malloc(sizeof(double)*total_particles);

        dbg_x = (double *)malloc(sizeof(double)*total_particles);
        dbg_y = (double *)malloc(sizeof(double)*total_particles);
        dbg_z = (double *)malloc(sizeof(double)*total_particles);

	 	readFortranBinaryArray(f,dbg_x);
	 	readFortranBinaryArray(f,dbg_y);
	 	readFortranBinaryArray(f,dbg_z);

	     for(int i = 0;i< total_particles;i++)
	     {

	          //px  = atof(str + 75);
	          //py  = atof(str + 100);
	          //pz  = atof(str + 125);
	          //m   = atof(str + 150);
	          //q_m = atof(str + 175);

		  Particle & p = vp[n++];
		//  p.SetXnext(x);
	     }
	}

	int PeriodicCurrentBoundaries(double* E,double *dbg_E, int dirE,int dir, int start1, int end1, int start2, int end2)
	{
	      Cell<Particle,dims>  c = (*AllCells)[0];

	      int N = getBoundaryLimit(dir);

	      for(int k = start2;k <= end2;k++)
	      {
	    	  for(int i = start1;i <= end1;i++)
		  {
		      int n1    = c.getGlobalBoundaryCellNumber(i,k,dir,1);
		      int n_Nm1 = c.getGlobalBoundaryCellNumber(i,k,dir,N-1);
	#ifdef DEBUG_PLASMA
		      int3 n1_3 = c.getCellTripletNumber(n1);
		      int3 n_Nm1_3 = c.getCellTripletNumber(n_Nm1);
	#endif
		      if(dir != dirE)
		      {
		         E[n1] += E[n_Nm1];
#ifdef PERIODIC_DEBUG_PRINTS
			      printf("to (%d,%d,%d) is added (%d,%d,%d): %15.5e += %15.5e \n",
			    		  n1_3.x+1,n1_3.y+1,n1_3.z+1,n_Nm1_3.x+1,n_Nm1_3.y+1,n_Nm1_3.z+1,E[n1],E[n_Nm1]);
#endif

		      }
		      if(dir != 1 || dirE != 1)
		      {
		         E[n_Nm1] =  E[n1];
		      }
	#ifdef DEBUG_PLASMA
	              int n1_check = CheckValue(E,dbg_E,n1);
		      if(n1_check == 0)
		      {
		      }
	#endif

		      int n_Nm2 = c.getGlobalBoundaryCellNumber(i,k,dir,N-2);
		      int n0    = c.getGlobalBoundaryCellNumber(i,k,dir,0);
	#ifdef DEBUG_PLASMA
		      int3 n0_3 = c.getCellTripletNumber(n0);
		      int3 n_Nm2_3 = c.getCellTripletNumber(n_Nm2);
	#endif

		      E[n0] += E[n_Nm2];
#ifdef PERIODIC_DEBUG_PRINTS
		      printf("to (%d,%d,%d) is added (%d,%d,%d): %15.5e += %15.5e \n",
		    		  n0_3.x+1,n0_3.y+1,n0_3.z+1,n_Nm2_3.x+1,n_Nm2_3.y+1,n_Nm2_3.z+1,E[n0],E[n_Nm2]);
#endif
		      E[n_Nm2] = E[n0];
	#ifdef DEBUG_PLASMA
	              int n_Nm2_check = CheckValue(E,dbg_E,n_Nm2);

	#endif

		      //   E[n0] = E[n_Nm2];
		      //   E[n_Nm1] = E[n1];

	#ifdef DEBUG_PLASMA
	              int n0_check    = CheckValue(E,dbg_E,n0);
	              int n_Nm1_check = CheckValue(E,dbg_E,n_Nm1);
	#endif

		     // }
		  }
	      }
	      return 0;
	}

	void ClearAllParticles(void )
	{
	    for(int n = 0;n < (*AllCells).size();n++)
	    {
	        Cell<Particle,dims> c = (*AllCells)[n];

		c.ClearParticles();

	    }
	}





	public:

	  GPUPlasma(int nx,int ny,int nz,double lx,double ly,double lz,
			  int n_per_cell1,double q_m,double TAU,int f3D,double B0,int N0,
			  double tex01,double tey01,double tez01,double Tb1,double rimp1,double rbd1,double ni1,
			  double bx,double by,double bz,
			  int beam_plasma1,int start_from_file)
	   {
		         beam_plasma = beam_plasma1;
		         nt_start_from_file = start_from_file;

		  	     tex0 = tex01;
		  	     tey0 = tey01;
		  	     tez0 = tez01;

		  	     Tb   = Tb1;
		  	     rimp = rimp1;
		  	     rbd  = rbd1;
                 ni = ni1;

		  	     lp   = n_per_cell1;

	     mesh.x = nx;
	     mesh.y = ny;
	     mesh.setZ(nz);

	     xmax.x = lx;
	     xmax.y = ly;
	     xmax.setZ(lz);

	     beam_max.x = bx;
	     beam_max.y = by;
	     beam_max.setZ(bz);

//	     ni = ni1;

	     n_per_cell = n_per_cell1;
	     ion_q_m    = q_m;
	     tau        = TAU;

	     flag3D  = f3D;

	     Bx0     = B0;
#ifdef DEBUG_WHOLE_PARTICLE_LIST
	     N = N0;
#else
	     N = mesh.size2()*n_per_cell;
#endif

	     beam_length = bx;
	     beam_width  = by;
	   }


	   virtual void InitializeCPU()
	   {
	      vector<Particle> vp;

	      Alloc();
	 //     exit(0);
	      Cell<Particle,dims> c000;

	      tgt_right  = new ParticleTarget<Particle,DIMENSIONS>;

	      tgt_right->Init(&(initial[2]),xmax.x,mesh);


	      InitCells();
	      c000 = (*AllCells)[0];
	//      exit(0);

	      InitFields();
	      c000 = (*AllCells)[0];
	      InitCurrents();

	      LoadTestData(START_STEP_NUMBER,START_STEP_NUMBER);
	      c000 = (*AllCells)[0];
	      magf = 1;

	      int size = mesh.size2();

	      cp = (Cell<Particle,dims> **)malloc(size*sizeof(Cell<Particle,dims> *));

	      for(int i = 0;i < size;i++)
	      	{
	      	 	Cell<Particle,dims> c,*d_c;
	      	// 	z0 = h_CellArray[i];
	      	    d_c = c.allocateCopyCellFromDevice();

	      	    cp[i] = d_c;
	      	}

	   }

	   void Free();

	   virtual void SetInitialConditions(){}

	   virtual void ParticleSort(){}

	   //void ApplyToAllParticles(cell_work_function cwf);

	//   virtual void ComputeField(int nt)
	//{
	//   double t27 = Hy[27];
	//
	//#ifdef UNITY_ELECTRIC_FIELD
	//     for(int i = 0;i < mesh.size2();i++)
	//     {
	//         Ex[i] = 1.0;
	//         Ey[i] = 1.0;
	//         Ez[i] = 1.0;
	//     }
	//#else
	//        double t271 = Hy[27];
	//
	//     CheckArray(Hx,dbgHx);
	//     CheckArray(Hy,dbgHy);
	//     CheckArray(Hz,dbgHz);
	//     CheckArray(Ex,dbgEx);
	//     CheckArray(Ey,dbgEy);
	//     CheckArray(Ez,dbgEz);
	//
	//     emh1(Qx,Qy,Qz,Hx,Hy,Hz,nt,Ex,Ey,Ez);
	//
	//     eme(Ex,Ey,Ez,nt,Hx,Hy,Hz,npJx,npJy,npJz);
	//     CheckArray(Ex,dbgEx);
	//     CheckArray(Ey,dbgEy);
	//     CheckArray(Ez,dbgEz);
	//
	//     emh2(Hx,Hy,Hz,nt+1,Qx,Qy,Qz);
	//
	//
	//     eme(Ex,Ey,Ez,nt+1,Hx,Hy,Hz,npJx,npJy,npJz);
	//
	//     CheckArray(Hx,dbgHx);
	//     CheckArray(Hy,dbgHy);
	//     CheckArray(Hz,dbgHz);
	//     puts("H");
	//     //exit(0);
	//
	////for(int i = 0;i < mesh.size2();i++)
	////     {
	////         Ex[i] = 0.0;
	////         Ey[i] = 0.0;
	////         Ez[i] = 0.0;
	////     }
	////    LoadTestData(3);
	//     eme(Ex,Ey,Ez,nt,Hx,Hy,Hz,npJx,npJy,npJz);
	//     CheckArray(Ex,dbgEx);
	//     CheckArray(Ey,dbgEy);
	//     CheckArray(Ez,dbgEz);
	//#endif
	//}

	//    void  ComputeField_FirstHalfStep(
	//		   double *locEx,double *locEy,double *locEz,
	//		   double nt,
	//		   double *locHx,double *locHy,double *locHz,
	//		   double *loc_npJx,double *loc_npJy,double *loc_npJz)
	//{
	//	 double *locQx,*locQy,*locQz;
	//	 static int first = 0;
	//
	//	 if(first == 0)
	//	 {
	//		 locQx = (double *)malloc(sizeof(double)*mesh.size2());
	//		 locQy = (double *)malloc(sizeof(double)*mesh.size2());
	//		 locQz = (double *)malloc(sizeof(double)*mesh.size2());
	//
	//		 first = 1;
	//	 }
	//     emh1(locQx,locQy,locQz,locHx,locHy,locHz,nt,locEx,locEy,locEz);
	//     CheckArray(dbg_Qx,locQx);
	//     eme(locEx,locEy,locEz,nt,locHx,locHy,locHz,loc_npJx,loc_npJy,loc_npJz);
	//     CheckArray(dbg_Qx,locQx);
	//     emh2(locHx,locHy,locHz,nt,locQx,locQy,locQz);
	//}




	   void ListAllParticles(int nt,string where)
	   	{
	   		FILE *f;
	   		char str[200];
	   		//Cell<Particle,dims> **cp;

//	   		printf("in ListAllParticles rank %5d \n",getRank());

	   		if(nt % minor_steps != 0)
	   		{
	   			return;
	   		}
	   		else
	   		{
#ifndef LIST_ALL_PARTICLES
	   		return;

#endif
	   		}

	   		sprintf(str,"List%05d_rank%05d_%s.dat",nt,getRank(),where.c_str());
//	   		puts(str);

	   		if((f = fopen(str,"wt")) == NULL) return;
//	   		puts(str);


	   		int size = (*AllCells).size();


	   		copyCells(where,nt);

	   			//h_ctrl = new Cell<Particle,dims>;

	   		for(int i = 0;i < size;i++)
	   		{
	   		 	Cell<Particle,dims> c = (*AllCells)[i];

   			    Particle p;
	 //   			int j;

	   			c.readParticleFromSurfaceDevice(0,&p);
	   	        //h_c.copyCellFromDevice(&d_c);
	   	        c.printFileCellParticles(f,cp[i]);
	   		}
	   		fclose(f);
	   	}

	void FortranOrder_StepAllCells(int nt)
	{
		int cell_sum = 0;
		int part_number = 0;
	//	double t_hx,t_hy,t_hz,*dbg;

		memset(Jx,0,sizeof(double)*mesh.size2());
		memset(Jy,0,sizeof(double)*mesh.size2());
		memset(Jz,0,sizeof(double)*mesh.size2());

		for(int n = 0;n < (*AllCells).size();n++)
		{
			Cell<Particle,dims> c = (*AllCells)[n];
	        part_number += c.number_of_particles;
		}

		for(int j = 1;j <= part_number;j++)
		{
		    for(int n = 0;n < (*AllCells).size();n++)
		    {
			    int f;

		        Cell<Particle,dims> c = (*AllCells)[n];

	            for(int i = 0; i < c.number_of_particles;i++)
	            {
	            	if(c.getFortranParticleNumber(i) != j) continue;
	            	int num = c.getFortranParticleNumber(i);
	            	c.SetAllCurrentsToZero();
		            f = c.Move(i);
		        //    sum += f;
		            c.writeAllToArrays(Jx,Jy,Jz,Rho,c.getFortranParticleNumber(i));
		            printf("particle number %10d \n",num);


	            }
		    }

	        //if(sum != c.number_of_particles)
		}
		printf("passed %10d cells of %10d total \n",cell_sum,(*AllCells).size());

		checkNonPeriodicCurrents(nt);

	    SetPeriodicCurrents(nt);
		CheckArray(Jx,dbgJx);
		CheckArray(Jy,dbgJy);
		CheckArray(Jz,dbgJz);

		//AssignCellsToArrays();
	}

	double TryCheckCurrent(int nt,double *npJx)
	{
		double *dbg,t_hx;//,t_hy,t_hz;



	  	dbg = (double *)malloc(sizeof(double)*mesh.size2());

	  	// read magnetic field from "nt+1" exlg file - to consider emh2



	    return 1.0;//t_hx;
	}

	double checkNonPeriodicCurrents(int nt)
	{

		printf("CHECKING Non-periodic currents !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n");

	  //	dbg = (double *)malloc(sizeof(double)*mesh.size2());
	  	TryCheckCurrent(nt,npJx);

	  	// read magnetic field from "nt+1" exlg file - to consider emh2

		return 1.0;//(t_hx+t_hy+t_hz)/3.0;
	}

	double checkPeriodicCurrents(int nt)
	{
		double *dbg,t_hx,t_hy,t_hz;

		printf("CHECKING periodic currents !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n");

	  	dbg = (double *)malloc(sizeof(double)*mesh.size2());

	  	TryCheckCurrent(nt,Jx);

	  	// read magnetic field from "nt+1" exlg file - to consider emh2

		return 1.0;//(t_hx+t_hy+t_hz)/3.0;
	}

	void printCurrentTensor(CellDouble dbg_cell_Jx,CellDouble dbg_cell_Jy,CellDouble dbg_cell_Jz,CurrentTensor t1)
	{
		double t_b,t_a,t;
		puts("Jx");
		t_b = dbg_cell_Jx.M[t1.Jx.i11][t1.Jx.i12][t1.Jx.i13];
		t   = t1.Jx.t[0];
		t_a = t_b + t;
		printf("before %15.5e t1 %15.5e after %15.5e \n",t_b,t,t_a);

		t_b = dbg_cell_Jx.M[t1.Jx.i21][t1.Jx.i22][t1.Jx.i23];
		t   = t1.Jx.t[1];
		t_a = t_b + t;
		printf("before %15.5e t1 %15.5e after %15.5e \n",t_b,t,t_a);

		t_b = dbg_cell_Jx.M[t1.Jx.i31][t1.Jx.i32][t1.Jx.i33];
		t   = t1.Jx.t[2];
		t_a = t_b + t;
		printf("before %15.5e t1 %15.5e after %15.5e \n",t_b,t,t_a);

		t_b = dbg_cell_Jx.M[t1.Jx.i41][t1.Jx.i42][t1.Jx.i43];
		t   = t1.Jx.t[3];
		t_a = t_b + t;
		printf("before %15.5e t1 %15.5e after %15.5e \n",t_b,t,t_a);

puts("Jy");
		printf("before %15.5e t1 %15.5e after %15.5Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½\n",
        dbg_cell_Jy.M[t1.Jy.i11][t1.Jy.i12][t1.Jy.i13],t1.Jy.t[0],dbg_cell_Jy.M[t1.Jy.i11][t1.Jy.i12][t1.Jy.i13]);
		printf("before %15.5e t1 %15.5e after %15.5Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½\n",
        dbg_cell_Jy.M[t1.Jy.i21][t1.Jy.i22][t1.Jy.i23],t1.Jy.t[1],dbg_cell_Jy.M[t1.Jy.i21][t1.Jy.i22][t1.Jy.i23]);
		printf("before %15.5e t1 %15.5e after %15.5Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½\n",
        dbg_cell_Jy.M[t1.Jy.i31][t1.Jy.i32][t1.Jy.i33],t1.Jy.t[2],dbg_cell_Jy.M[t1.Jy.i31][t1.Jy.i32][t1.Jy.i33]);
		printf("before %15.5e t1 %15.5e after %15.5Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½\n",
        dbg_cell_Jy.M[t1.Jy.i41][t1.Jy.i42][t1.Jy.i43],t1.Jy.t[3],dbg_cell_Jy.M[t1.Jy.i41][t1.Jy.i42][t1.Jy.i43]);
        puts("Jz");
		printf("before %15.5e t1 %15.5e after %15.5Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½\n",
        dbg_cell_Jz.M[t1.Jz.i11][t1.Jz.i12][t1.Jz.i13],t1.Jz.t[0],dbg_cell_Jz.M[t1.Jz.i11][t1.Jz.i12][t1.Jz.i13]);
		printf("before %15.5e t1 %15.5e after %15.5Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½\n",
        dbg_cell_Jz.M[t1.Jz.i21][t1.Jz.i22][t1.Jz.i23],t1.Jz.t[1],dbg_cell_Jz.M[t1.Jz.i21][t1.Jz.i22][t1.Jz.i23]);
		printf("before %15.5e t1 %15.5e after %15.5Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½\n",
        dbg_cell_Jz.M[t1.Jz.i31][t1.Jz.i32][t1.Jz.i33],t1.Jz.t[2],dbg_cell_Jz.M[t1.Jz.i31][t1.Jz.i32][t1.Jz.i33]);
		printf("before %15.5e t1 %15.5e after %15.5Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½Ã¯Â¿Â½\n",
        dbg_cell_Jz.M[t1.Jz.i41][t1.Jz.i42][t1.Jz.i43],t1.Jz.t[3],dbg_cell_Jz.M[t1.Jz.i41][t1.Jz.i42][t1.Jz.i43]);
	}


	int write_field_component(int nt,double *d_f,string where,string name,int size)
	{
		static int first = 1;
		static double *h_f;

#ifndef DEBUG_FIELD_COMPONENT_PRINTS
		return 0;
#endif

		if(first == 1)
		{
			first = 0;

			MemoryAllocate((void **)&h_f,sizeof(double)*size);
		}

		MemoryCopy(h_f,d_f,sizeof(double)*size,DEVICE_TO_HOST);

		 Cell<Particle,dims> c = (*AllCells)[0];
		 double hx = c.get_hx(),hy = c.get_hy(),hz = c.get_hz();
		 char fname[100];
		 FILE *f;

		 sprintf(fname,"comp%s_at_%s_%s_rank%05d_nt%08d.dat",name.c_str(),where.c_str(),unique_variant_name,getRank(),nt);

		 if((f = fopen(fname,"wt")) == NULL) return 1;

		 for(int i = 0;i < c.mesh.x;i++)
		 {
			 for(int l = 0;l < c.mesh.y;l++)
			 {
				 for(int k = 0;k < c.mesh.z();k++)
				 {
					 int n  = c.getGlobalCellNumber(i,l,k);

					 fprintf(f,"%10d %5d %5d %10d %10.3e %10.3e %10.3e %25.15e \n",i,l,k,n,i*hx,l*hy,k*hz,h_f[n]);
				 }
			 }
		 }
		 fclose(f);


		 return 0;
	}



	void CellOrder_StepAllCells(int nt,double mass,double q_mass,int first)
	{

//		dim3 dimGrid(mesh.x+2,mesh.y+2,mesh.dimz2()),dimGridOne(1,1,1),dimBlock(512,1,1),
//				dimBlockOne(1,1,1),dimBlockGrow(1,1,1),dimBlockExt(CellExtent,CellExtent,CellExtent);
//		dim3 dimGridBulk(mesh.x,mesh.y,mesh.z());
		int  err;

		memory_monitor("CellOrder_StepAllCells1",nt);
//		write_field_component(nt,d_Jx,"beginCellOrder","Jx",mesh.size2());


		memset(Jx,0,sizeof(double)*mesh.size2());
		memset(Jy,0,sizeof(double)*mesh.size2());
		memset(Jz,0,sizeof(double)*mesh.size2());
		MemorySet(d_Jx,0,sizeof(double)*mesh.size2());
		MemorySet(d_Jy,0,sizeof(double)*mesh.size2());
	 	MemorySet(d_Jz,0,sizeof(double)*mesh.size2());


		char name[100];

		sprintf(name,"before_set_to_zero_%03d.dat",nt);
		write3D_GPUArray(name,d_Jx);


    	    SetAllCurrentsToZero(d_CellArray,mesh.x+2,mesh.y+2,mesh.dimz2(),
    	    		CellExtent,CellExtent,CellExtent);

//	        GPU_SetAllCurrentsToZero<<<dimGrid, dimBlockExt,16000>>>(d_CellArray);
		    err = getLastError();
//		    ThreadSynchronize();
//		    TEST_ERROR(err);
		    memory_monitor("CellOrder_StepAllCells3",nt);

			sprintf(name,"before_step_%03d.dat",nt);
			write3D_GPUArray(name,d_Jx);
//			printf("b ListAllParticles rank %5d \n",getRank());
			ListAllParticles(nt,"bStepAllCells");

			params.mass = mass;
			params.q_mass = q_mass;
			params.jmp    = jmp;
			params.nt     = nt;
			params.d_ctrlParticles = d_ctrlParticles;
			params.particles_processed_by_a_single_thread = PARTICLE_BLOCK_X_SIZE;
			params.flown_beam_particles  = 0;


			MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

		    Kernel_Launcher(d_CellArray,d_params,mesh.x+2,mesh.y+2,mesh.dimz2(),
		    		params.particles_processed_by_a_single_thread,
		    		1,1,16000,h_GPU_StepAllCells_SingleNode,
//		    		PARTICLE_BLOCK_X_SIZE,1,1,16000,h_GPU_StepAllCells_SingleNode,
		            "StepAllCells");

//		    write_field_component(nt,d_Jx,"afterStepAllCells","Jx",mesh.size2());

		    DeviceSynchronize();

		    MemoryCopy(&params,d_params,sizeof(KernelParams),DEVICE_TO_HOST);

		    flown_beam_particles = params.flown_beam_particles;

//            GPU_StepAllCells<<<dimGrid, dimBlock,16000>>>(d_CellArray/*,d_jx,d_jy,d_jz*/,
//            		     		                          mass,q_mass,d_ctrlParticles,jmp,nt);
            
            err = getLastError();
            ThreadSynchronize();
//            TEST_ERROR(err);

            memory_monitor("CellOrder_StepAllCells4",nt);

            //ListAllParticles(nt,"aStepAllCells");


		    int err1 = getLastError();
		    DeviceSynchronize();
		    int err2 = getLastError();
		    char err_s[200];
		    strcpy(err_s,getErrorString(err2));


			sprintf(name,"before_write_currents_%03d.dat",nt);
			write3D_GPUArray(name,d_Jx);

//		    dim3 dimExt(CellExtent,CellExtent,CellExtent);


		    params.d_Jx = d_Jx;
	        params.d_Jy = d_Jy;
		   	params.d_Jz = d_Jz;

		   	params.d_Rho = d_Rho;
		    	//		static int first  = 1;
	    MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

//	    write_field_component(nt,d_Jx,"bWrite","Jx",mesh.size2());
    	Kernel_Launcher(d_CellArray,d_params,mesh.x+2,mesh.y+2,mesh.dimz2(),
    			CellExtent,CellExtent,CellExtent,16000,
    	                h_GPU_WriteAllCurrents_SingleNode,
    	               "WriteAllCurrents");
//    	write_field_component(nt,d_Jx,"awrite","Jx",mesh.size2());


//	    GPU_WriteAllCurrents<<<dimGrid, dimExt,16000>>>(d_CellArray,0,d_Jx,d_Jy,d_Jz,d_Rho);
 		    err = getLastError();
 		    ThreadSynchronize();
// 		    TEST_ERROR(err);

 		   memory_monitor("CellOrder_StepAllCells5",nt);

 			            sprintf(name,"after_write_currents_%03d.dat",nt);
 						write3D_GPUArray(name,d_Jx);
#ifdef PRINT_CELL_CURRENTS
#endif
// 						write_field_component(nt,d_Jx,"bSUM","Jx",mesh.size2());

 						MemoryCopy(Jx,d_Jx,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
 						MemoryCopy(Jy,d_Jy,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
 						MemoryCopy(Jz,d_Jz,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);

 						sumMPI(mesh.size2(),Jx,Jy,Jz);

 						MemoryCopy(d_Jx,Jx,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
 						MemoryCopy(d_Jy,Jy,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);
 						MemoryCopy(d_Jz,Jz,sizeof(double)*mesh.size2(),HOST_TO_DEVICE);

 						write_field_component(nt,d_Jx,"aSUM","Jx",mesh.size2());
 						write_field_component(nt,d_Ex,"aSUM","Ex",mesh.size2());


 						memory_monitor("CellOrder_StepAllCells6",nt);

 						int before_MakeDepartureLists,after_MakeDepartureLists,
                          before_ArrangeFlights,after_ArrangeFlights;


#ifdef BALANCING_PRINTS
              before_MakeDepartureLists = getLastError();
              printf("before_MakeDepartureLists %d %s blockdim %d %d %d\n",before_MakeDepartureLists,
            		  GetErrorString(before_MakeDepartureLists),dimGrid.x,dimGrid.y,dimGrid.z);
#endif

              int *stage,*stage1,*d_stage,*d_stage1;

              stage  = (int *)malloc(sizeof(int)*mesh.size2());
              stage1 = (int *)malloc(sizeof(int)*mesh.size2());

              MemoryAllocate((void**)&d_stage,sizeof(int)*mesh.size2());
              MemoryAllocate((void**)&d_stage1,sizeof(int)*mesh.size2());

          	  params.nt = nt;
          	  params.d_stage = d_stage;
          	  MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

          	  Kernel_Launcher(d_CellArray,d_params,mesh.x+2,mesh.y+2,mesh.dimz2(),
          	                1,1,1,16000,h_GPU_MakeDepartureLists_SingleNode,
          	                "MakeDepartureLists");

//              GPU_MakeDepartureLists<<<dimGrid, dimBlockOne>>>(d_CellArray,nt,d_stage);
              after_MakeDepartureLists = getLastError();
              err = getLastError();
              ThreadSynchronize();
//              TEST_ERROR(after_MakeDepartureLists);
#ifdef BALANCING_PRINTS
              printf("after_MakeDepartureLists %d %s\n",after_MakeDepartureLists,GetErrorString(after_MakeDepartureLists));
#endif

                                DeviceSynchronize();
              err = MemoryCopy(stage,d_stage,sizeof(int)*mesh.size2(),DEVICE_TO_HOST);
              if(err != 0)
              {
            	  puts("copy error");
            	  exit(0);
              }
              if(stage[0] == TOO_MANY_PARTICLES)
              {
            	  printf("too many particles flying to (%d,%d,%d) from (%d,%d,%d) \n",
            			  stage[1],stage[2],stage[3],
            			  stage[4],stage[5],stage[6]);
            	  exit(0);
              }

             // ListAllParticles(nt,"aMakeDepartureLists");
              //exit(0);
#ifdef BALANCING_PRINTS
              before_ArrangeFlights = getLastError();
              printf("before_ArrangeFlights %d %s\n",before_ArrangeFlights,GetErrorString(before_ArrangeFlights));
#endif


              MemorySet(d_stage1,0,sizeof(int)*mesh.size2());
          	  params.d_stage = d_stage1;
          	  MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

          	  Kernel_Launcher(d_CellArray,d_params,mesh.x+2,mesh.y+2,mesh.dimz2(),
          	                1,1,1,16000,h_ArrangeFlights,
          	                "ArrangeFlights");


//              GPU_ArrangeFlights<<<dimGridBulk, dimBlockOne>>>(d_CellArray,nt,d_stage1);
              after_ArrangeFlights = getLastError();
              ThreadSynchronize();
//              TEST_ERROR(after_MakeDepartureLists);
#ifdef BALANCING_PRINTS
              printf("after_ArrangeFlights %d %s\n",after_ArrangeFlights,GetErrorString(after_ArrangeFlights));
                                DeviceSynchronize();
#endif

              err = MemoryCopy(stage1,d_stage1,sizeof(int)*mesh.size2(),DEVICE_TO_HOST);
                                              if(err != 0)
                                              {
                                            	  puts("copy error");
                                            	  exit(0);
                                              }
             // ListAllParticles(nt,"aArrangeFlights");


              memory_monitor("CellOrder_StepAllCells7",nt);


	}

	void SetCurrentsFromCellsToArrays(int nt)
	{
	    memset(Jx,0,sizeof(double)*mesh.size2());
	    memset(Jy,0,sizeof(double)*mesh.size2());
	    memset(Jz,0,sizeof(double)*mesh.size2());

	    for(int n = 0;n < (*AllCells).size();n++)
	    {

	        Cell<Particle,dims> c = (*AllCells)[n];

	        c.writeAllToArrays(Jx,Jy,Jz,Rho,c.getFortranParticleNumber(0));

	    }

	    //if(sum != c.number_of_particles)
	    //	}

	    memcpy(npJx,Jx,sizeof(double)*mesh.size2());
	   	memcpy(npJy,Jy,sizeof(double)*mesh.size2());
	   	memcpy(npJz,Jz,sizeof(double)*mesh.size2());
	   	checkNonPeriodicCurrents(nt);
	    SetPeriodicCurrents(nt);
	    CheckArray(Jx,dbgJx);
	    CheckArray(Jy,dbgJy);
	    CheckArray(Jz,dbgJz);

	}

	//   virtual void Diagnose(){ puts("Plasma");}

	   void virtual Step()
	{
	    // double t = dbgJx[0];
	  //   int j;
	     Particle p;

	     //ComputeField();
		   Cell<Particle,dims> & c000 = (*AllCells)[0];
	    // AssignArraysToCells();

	     StepAllCells();

	     puts("particles moved!!!");
	     //return;
	     exit(0);
	/*
	     c000 = (*AllCells)[0];
	     AssignCellsToArrays();

	     for(int n = 0;n < (*AllCells).size();n++)
	     {
	         Cell<Particle,dims>  c = (*AllCells)[n];
		     thrust::host_vector<Particle> vecp = c.getFlyList();
		     if(vecp.size() > 0)
		     {
		        int q = 0;
		        for(j = 0;j < vecp.size();j++)
		        {
		        	p = vecp[j];
		        	printf("cell %5d particle %5d %10.3e \n",n,j,p.x);
		        }
		     }
		     Distribute(vecp);
	     }
	     c000 = (*AllCells)[0];
	   //  LoadTestData(2);

	     ComputeField();

	#ifdef DEBUG_PLASMA
	     CheckArray(Jx, dbgJx);
	#endif

	     ParticleLog();*/
	}

//double checkControlPointParticlesOneSort(int check_point_num,FILE *f,GPUCell<Particle,dims> **copy_cells,int nt,int sort)
//{
//
//    double t = 0.0;
//    int size = 1;
//#ifdef CPU_DEBUG_RUN
//    double q_m,m;
//    struct sysinfo info;
//
//    memory_monitor("checkControlPointParticlesOneSort",nt);
//
//  //  double x,y,z,px,pz,q_m,*buf,tp,m;
//    //double dbg_x,dbg_y,dbg_z,dbg_px,dbg_py,dbg_pz;
//
//    Cell<Particle,dims> c0 = (*AllCells)[0];
//    //int pn_min/*,pn_ave,pn_max*/;
//    if(check_point_num == 50)
//    {
//    	int qq = 0;
//    }
//    total_particles = readBinaryParticleArraysOneSort(f,&dbg_x,&dbg_y,&dbg_z,
//   		                                             &dbg_px,&dbg_py,&dbg_pz,&q_m,&m,nt,sort);
//    memory_monitor("checkControlPointParticlesOneSort2",nt);
//
//    size = (*AllCells).size();
//
//   	for(int i = 0;i < size;i++)
//   	{
//   	 	GPUCell<Particle,dims> c = *(copy_cells[i]);
//
//#ifdef checkControlPointParticles_PRINT
//             printf("cell %d particles %20d \n",i,c.number_of_particles);
//#endif
//
//   	 	t += c.checkCellParticles(check_point_num,dbg_x,dbg_y,dbg_z,dbg_px,dbg_py,dbg_pz,q_m,m);
//   	}
//   	memory_monitor("checkControlPointParticlesOneSort3",nt);
//
//   	free(dbg_x);
//   	free(dbg_y);
//   	free(dbg_z);
//
//   	free(dbg_px);
//   	free(dbg_py);
//   	free(dbg_pz);
//   	memory_monitor("checkControlPointParticlesOneSort4",nt);
//#endif
//	return t/size;
//}

double checkControlPointParticles(int check_point_num,FILE *f,char *fname,int nt)
{
	double te = 0.0,ti = 0.0,tb = 0.0;
	struct sysinfo info;
#ifdef CPU_DEBUG_RUN
 //   Cell<Particle,dims> **cp;

	int size = (*AllCells).size();


	char where[100];
	sprintf(where, "checkpoint%03d",check_point_num);
	copyCells(where,nt);

	//checkParticleNumbers(cp,check_point_num);

#ifdef FREE_RAM_MONITOR
	sysinfo(&info);
#ifdef checkControlPointParticles_PRINTS
	printf("checkControlPointParticles %u \n",info.freeram/1024/1024);
#endif
#endif


	GPUCell<Particle,dims> c = *(cp[141]);
#ifdef checkControlPointParticles_PRINTS
	printf("checkControlPointParticlesOneSort cell 141 particles %20d \n",c.number_of_particles);
#endif

#ifdef FREE_RAM_MONITOR
	sysinfo(&info);
//	printf("checkControlPointParticles0.9 %u \n",info.freeram/1024/1024);
#endif

	ti  = checkControlPointParticlesOneSort(check_point_num,f,cp,nt,ION);
#ifdef FREE_RAM_MONITOR
	sysinfo(&info);
//	printf("checkControlPointParticles1 %u \n",info.freeram/1024/1024);
#endif

	te  = checkControlPointParticlesOneSort(check_point_num,f,cp,nt,PLASMA_ELECTRON);

#ifdef FREE_RAM_MONITOR
	sysinfo(&info);
//    printf("checkControlPointParticles1.5 %u \n",info.freeram/1024/1024);
#endif

	tb  = checkControlPointParticlesOneSort(check_point_num,f,cp,nt,BEAM_ELECTRON);

#ifdef FREE_RAM_MONITOR
	sysinfo(&info);
//	printf("checkControlPointParticles2 %u \n",info.freeram/1024/1024);
#endif




#endif

    memory_monitor("after_free",nt);
	return (te+ti+tb)/3.0;
}

int readControlFile(int nt)
{


#ifndef ATTRIBUTES_CHECK
	return 0;
#else
	FILE *f;
	char fname[100];
	static int first = 1;
	int size;//,jmp1;

	sprintf(fname,"ctrl%05d",nt);

	if((f = fopen(fname,"rb")) == NULL)
		{
		  puts("no ini-file");
		  exit(0);
		}

	fread(&size,sizeof(int),1,f);
	fread(&ami,sizeof(double),1,f);
	fread(&amf,sizeof(double),1,f);
	fread(&amb,sizeof(double),1,f);
	fread(&size,sizeof(int),1,f);

	fread(&size,sizeof(int),1,f);

	if(first == 1)
	{
		first = 0;
        ctrlParticles = (double *)malloc(size);
#ifdef ATTRIBUTES_CHECK
        memset(ctrlParticles,0,size);
        MemoryAllocate((void**)&d_ctrlParticles,size);
        MemorySet(d_ctrlParticles,0,size);
        size_ctrlParticles = size;
#endif
	}
	fread(ctrlParticles,1,size,f);


	jmp = size/sizeof(double)/PARTICLE_ATTRIBUTES/3;

	//double x,y,z;
	//int pos;

	//pos = ParticleAttributePositionFortran(jmp,1,ION,1);
	//x = ctrlParticles[pos];
	//printf("x %s \n",FortranExpWrite(x));

	//pos = ParticleAttributePositionFortran(jmp,1,ION,2);

	//y = ctrlParticles[pos];
	//pos = ParticleAttributePositionFortran(jmp,1,ION,3);

	//z = ctrlParticles[pos];

	return 0;
#endif
}

int checkDoublePrecisionIdentity(double a,double b)
{
	char as[50],bs[50];
	int point_pos,i;

	if(fabs(a-b) > PARTICLE_TOLERANCE) return 0;

	sprintf(as,"%25.15e",a);
	sprintf(bs,"%25.15e",b);

	for(i = 0;i < strlen(as);i++)
	{
		if(as[i] == bs[i] && as[i] == '.')
		{
			point_pos = i;
		}

		if(as[i] != bs[i]) break;
	}


	return (i - point_pos);
}

double checkParticleSortAttributes(int nt,particle_sorts sort,int attributes_checked,int jmp_real)
{


#ifndef ATTRIBUTES_CHECK
	return 0;
#else
	int j,i,n,n_fortran/*,n1,n1_fortran,n2,n2_fortran,max_j,max_num*/,eq_flag,wrong = 0;
		int min_eq_flag = 25;
		double t = 0.0,c,ch, delta,max_delta = 0.0;
		//c1,ch1,c2,ch2,
		double wrong_attr;
		FILE *f,*f_all,*f_tab;
		char fn[200],fn_all[200],fn_table[200];

	sprintf(fn,"attributes%05d_%d.dat",nt,(int)sort);
	sprintf(fn_all,"atr_all_attribute%05d_%d.dat",nt,(int)sort);
	sprintf(fn_table,"atr_table_attribute%05d_%d.dat",nt,(int)sort);

	if((f = fopen(fn,"wt")) == NULL) return -1.0;
	if((f_all = fopen(fn_all,"wt")) == NULL) return -1.0;
	if((f_tab = fopen(fn_table,"wt")) == NULL) return -1.0;

	for(i = 1;
			i <= attributes_checked
	//PARTICLE_ATTRIBUTES
	;i++)
	{
		wrong_attr = 0.0;
		max_delta  = 0.0;
		min_eq_flag = 25;
		for(j = 1;j <= jmp_real;j++)
		{
			n         = ParticleAttributePosition(jmp,j,sort,i);
			n_fortran = ParticleAttributePositionFortran(jmp,j,sort,i);


			c  = ctrlParticles[n_fortran];
			ch = check_ctrlParticles[n];
			delta = fabs(ctrlParticles[n_fortran] - check_ctrlParticles[n]);

			eq_flag = checkDoublePrecisionIdentity(c,ch);

			if(delta > max_delta)
				{
				max_delta = delta;
				}

			if(eq_flag < min_eq_flag) min_eq_flag = eq_flag;


			if(eq_flag >= TOLERANCE_DIGITS_AFTER_POINT) t += 1.0;
			else
			{
				wrong_attr += 1.0;
				eq_flag = checkDoublePrecisionIdentity(c,ch);
				fprintf(f,"wrong %10d particle %10d attribute %3d digits %5d CPU %25.16e GPU %25.16e diff %15.5e n %10d nf %10d \n",
						wrong++,j,i,eq_flag,c,ch,delta,n,n_fortran);
			}


		    fprintf(f_all,"wrong %10d particle %10d attribute %3d digits %5d CPU %25.16e GPU %25.16e diff %15.5e n %10d nf %10d \n",
				wrong++,j,i,eq_flag,c,ch,delta,n,n_fortran);

		}
		fprintf(f_tab,"sort %2d attribute %3d wrong %e, %8d of %10d delta %15.5e digits %5d\n",
				   (int)sort,i,wrong_attr/jmp_real,(int)wrong_attr,jmp_real,max_delta,min_eq_flag);
	}
	fclose(f);
	fclose(f_all);
	fclose(f_tab);
	return (1.0 + t/jmp/attributes_checked);
#endif
}

int checkParticleAttributes(int nt)
{


#ifndef ATTRIBUTES_CHECK
	return 0;
#else

	static int first = 1;

	readControlFile(nt);

	if(first == 1)
	{
		first = 0;
		check_ctrlParticles = (double *)malloc(size_ctrlParticles);
		memset(check_ctrlParticles,0,size_ctrlParticles);

	}
	int err;

	err = getLastError();

#ifdef ATTRIBUTES_CHECK
    err = MemoryCopy(check_ctrlParticles,d_ctrlParticles,
    		   //1,
    		   size_ctrlParticles, // /PARTICLE_ARRAY_PORTION,
    		   DEVICE_T_HOST);
    err = getLastError();
#endif
    if(err != 0)
    {
    	printf("MemoryCopy before attributes error %d %s\n",err,GetErrorString(err));
    	exit(0);
    }

    checkParticleSortAttributes(nt,ION,131,real_number_of_particle[(int)ION]);
    checkParticleSortAttributes(nt,PLASMA_ELECTRON,131,real_number_of_particle[(int)PLASMA_ELECTRON]);
    checkParticleSortAttributes(nt,BEAM_ELECTRON,131,real_number_of_particle[(int)BEAM_ELECTRON]);

#endif
}

void printGPUParticle(int num,int sort)
{
	int err;
//	dim3 dimGrid(mesh.x+2,mesh.y+2,mesh.dimz2()),dimGridOne(1,1,1),dimBlock(MAX_particles_per_cell/2,1,1),
//					dimBlockOne(1,1,1),dimBlockGrow(1,1,1),dimBlockExt(CellExtent,CellExtent,CellExtent);


//	printParticle<<<dimGrid, dimBlock,16000>>>(d_CellArray,num,sort);
    err = getLastError();
    ThreadSynchronize();
//    TEST_ERROR(err);
}





int checkParticleNumbers(GPUCell<Particle,dims> ** h_cp,int num)
{
	int *d_numbers,*h_numbers,size;
	int err;

	if(num >= 270) return -1;

	size = (*AllCells).size();

	h_numbers = (int *)malloc(size*sizeof(int));
	MemoryAllocate((void**)&d_numbers,size*sizeof(int));

	params.numbers = d_numbers;
	MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

	Kernel_Launcher(d_CellArray,d_params,mesh.size2(),1,1,
	                1,1,1,16000,h_GPU_GetCellNumbers_SingleNode,
	                "GetCellNumbers");

//	GPU_GetCellNumbers<<<mesh.size2(),1>>>(d_CellArray,d_numbers);
    err = getLastError();
    ThreadSynchronize();
//    TEST_ERROR(err);

	err = MemoryCopy(h_numbers,d_numbers,size*sizeof(int),DEVICE_TO_HOST);


	for(int i = 0;i < (*AllCells).size();i++)
	{
		GPUCell<Particle,dims> c = *(h_cp[i]);

	    if(h_numbers[i] != h_controlParticleNumberArray[i])
	    {
	    	printf("checkpoint %d: cell %d incorrect number of particles in DEVICE array %15d (must be %15d)\n",
	    			num,i,
	    			h_numbers[i],h_controlParticleNumberArray[i]
	    			);
	    	exit(0);
	    }
	}


	for(int i = 0;i < (*AllCells).size();i++)
	{
		GPUCell<Particle,dims> c = *(h_cp[i]);

	    if(c.number_of_particles != h_controlParticleNumberArray[i])
	    {
	    	printf("checkpoint %d: cell %d incorrect number of particles in HOST copy array %15d (must be %15d)\n",
	    			num,i,
	    			c.number_of_particles,h_controlParticleNumberArray[i]
	    			);
	    	exit(0);
	    }
	}

    return 0;
}

void printCellCurrents(int num,int nt,char *name,char *where)
{
	int size = (*AllCells).size();
	GPUCell<Particle,dims> ** cp = (GPUCell<Particle,dims> **)malloc(size*sizeof(GPUCell<Particle,dims> *));
	FILE *f;
	char fname[100];
	CellDouble *m;

	sprintf(fname,"%s_at_%s_cells_%03d.dat",name,where,nt);

	if((f = fopen(fname,"wt")) == NULL) return;

	//	copyCells(h_CellArray);
		copyCells(where,nt);

		for(int i = 0;i <size;i++)
		{
			GPUCell<Particle,dims> c = *(cp[i]);

			for(int i1  = 0;i1 < CellExtent;i1++)
			{
				for(int k1  = 0;k1 < CellExtent;k1++)
				{
					for(int l1  = 0;l1 < CellExtent;l1++)
					{
						if(!strcmp(name,"jx"))
						{
							m = c.Jx;
						}
						else
						{
							if(!strcmp(name,"jy"))
							{
								m = c.Jy;
							}
							else
							{
								m = c.Jz;
							}
						}
						fprintf(f,"%10d %5d %5d %5d %25.15e \n",i,i1,k1,l1,m->M[i1][k1][l1]);
					}
				}
			}
		}
        fclose(f);
}

int memory_monitor(string legend,int nt)
{
	static int first = 1;
	static FILE *f;

#ifndef FREE_RAM_MONITOR
	return 1;
#endif

	if(first == 1)
	{
		first = 0;
		f = fopen("memory_monitor.log","wt");
	}

	size_t m_free,m_total;
	struct sysinfo info;


	int err = GetDeviceMemory(&m_free,&m_total);

	sysinfo(&info);
	fprintf(f,"step %10d %50s GPU memory total %10lu free %10lu free CPU memory %10lu \n",nt,legend.c_str(),m_total/1024/1024,m_free/1024/1024,info.freeram/1024/1024);

}

void BeamInput(int nt)
{
     double dN = N;
     int d_N;
     double *x,*y,*z,*pu,*pv,*pw;
     Particle *p;

     dN /= mesh.size2();
     dN *= mesh.y*mesh.z();
     d_N = dN;

     KernelParams params;
//     params.flown_beam_particles  = 0;

//     MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);
//
//          Kernel_Launcher(d_CellArray,d_params,mesh.x,mesh.y+2,mesh.dimz2(),
//                  	                1,1,1,16000,h_GPU_GetFlownBeamNumber_SingleNode,
//                  	                "GPU_GetFlownBeamNumber_SingleNode");
//    MemoryCopy(&params,d_params,sizeof(KernelParams),DEVICE_TO_HOST);

     d_N = flown_beam_particles;
#ifdef BEAM_ADD_PRINTS
     printf("TOTAL FLOWN %d\n",d_N);
#endif



     vector<Particle> vp;
     double beam_q_m,beam_m;
     FILE *f;

     FILE *fa;
     char fname[100];

     sprintf(fname,"add_beam_stat_rank%05d_nt%08d.dat",getRank(),nt);


     beam_q_m = initial[2].q_m;
     beam_m   = initial[2].m[0];

     int jmb_real;
     double beam_hx = xmax.x/mesh.x;

     sumMPIint(&d_N);

      x  = (double *)malloc(d_N*sizeof(double));
      y  = (double *)malloc(d_N*sizeof(double));
      z  = (double *)malloc(d_N*sizeof(double));
      pu = (double *)malloc(d_N*sizeof(double));
      pv = (double *)malloc(d_N*sizeof(double));
      pw = (double *)malloc(d_N*sizeof(double));

     if((fa = fopen(fname,"wt")) == NULL) return;



     AddBeamParticles(d_N,
     				   tex0,tey0,tez0,
     				   beam_length, beam_width,&jmb_real,
                       xmax.x,xmax.y,xmax.z(),meh,Tb,rimp,rbd,
     				   x,y, z, pu,pv,pw);

     Particle *p_send_array,*d_p_send_array;
     p_send_array = (Particle *)malloc(sizeof(Particle)*d_N);

     int d_N_global = d_N;
     d_N = 0;

     for(int n = 0; n < d_N_global;n++)
     {
    	 p = new Particle(x[n],y[n],z[n],pu[n],pv[n],pw[n],beam_m,beam_q_m);
    	 p->sort = BEAM_ELECTRON;
    	 p->fortran_number = initial[2].total + n + 1;

    	 vp.push_back(*p);

       	 fprintf(fa,"total  %10d %5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",
        	          	    p->fortran_number,n,p->X.x,p->X.y,p->X.z(),p->pu,p->pv,p->pw);
     }


//     if((fa = fopen(fname,"at")) == NULL) return;

//     printf("d_N6 %10d rank %d \n",d_N,getRank());

     for(int n = 0; n < d_N_global;n++)
     {
//         printf("d_N7 %10d rank %d \n",d_N,getRank());

    	 *p = vp[n];
    	 if((p->fortran_number % getSize()) == getRank())
    	 {

    	     p_send_array[d_N] = *p;
    	     d_N++;
        	     fprintf(fa,"select %10d %5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",
        	          	        	p->fortran_number,n,
        	          	        		  p->X.x,p->X.y,p->X.z(),p->pu,p->pv,p->pw
        	          	        		  );


    	 }
//         printf("d_N8 %10d rank %d \n",d_N,getRank());

     }



//     printf("final d_N %10d rank %d \n",d_N,getRank());


//     if((fa = fopen(fname,"at")) == NULL) return;

     for(int n = 0; n < d_N;n++)
     {


         		      Particle *p = &(p_send_array[n]);

         	          fprintf(fa,"final  %10d %5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",
         	        		             p->fortran_number,n,
         	        		             p->X.x,p->X.y,p->X.z(),p->pu,p->pv,p->pw
         	        		  );


     }

     fclose(fa);

     MemoryAllocate((void **)&d_p_send_array,sizeof(Particle)*d_N);
     MemoryCopy(d_p_send_array,p_send_array,sizeof(Particle)*d_N,HOST_TO_DEVICE);



     params.d_p_send_array = (int *)d_p_send_array;
     params.p_send_array_size = d_N;
     MemoryCopy(d_params,&params,sizeof(KernelParams),HOST_TO_DEVICE);

     Kernel_Launcher(d_CellArray,d_params,mesh.x+2,mesh.y+2,mesh.dimz2(),
             	                1,1,1,16000,h_AddParticlesToCells,
             	                "AddParticlesToCells");



     initial[2].total += d_N_global;
     diagnostics[2].total += d_N_global;

}

};







#endif /* GPU_PLASMA_H_ */
