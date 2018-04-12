/*
 * gpucell.h
 *
 *  Created on: Oct 19, 2013
 *      Author: snytav
 */

#ifndef GPUCELL_H_
#define GPUCELL_H_

#include "cell.h"
#include "archAPI.h"



void dbgPrintGPUParticleAttribute(Cell<Particle,DIMENSIONS> *d_c,int n_particle,int attribute,char *name )
{
    double t;
    Cell<Particle,DIMENSIONS> *h_c;
    int shift = (attribute + n_particle*sizeof(Particle)/sizeof(double));
    int err;

    h_c = new Cell<Particle,DIMENSIONS>;

    err = MemoryCopy(h_c,d_c,sizeof(Cell<Particle,DIMENSIONS>),DEVICE_TO_HOST);
    if(err != 0)
        {
        	printf("pyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }
    double *vec = h_c->doubParticleArray + shift;

    MemoryCopy(&t,vec,sizeof(double),DEVICE_TO_HOST);

    printf("%s %10.3e \n",name,t);
}

//__global__ void testKernelBefore(double *vec,int n_particle,int attribute)
//{
//   //  	cuPrintf("vecBefore %15.5e \n",vec[attribute + n_particle*sizeof(Particle)/sizeof(double)]);
//}
//
//__global__ void testKernel(double *vec)
//{
//    // 	cuPrintf("vec %15.5e \n",vec[1]);
//}

template <class Particle,int dims >
class GPUCell: public Cell<Particle,dims>
{


public:
	  double *d_wrong_current_particle_attributes,*h_wrong_current_particle_attributes;



 #ifdef __CUDACC__
 __host__ __device__
 #endif

    GPUCell(){}

 #ifdef __CUDACC__
 __host__ __device__
 #endif

   ~GPUCell(){}

 #ifdef __CUDACC__
 __host__ __device__
 #endif

    GPUCell(int i1,int l1,int k1,double Lx,double Ly, double Lz,int Nx1, int Ny1, int Nz1,double tau1):
       Cell<Particle,dims>(i1,l1,k1,Lx,Ly,Lz,Nx1,Ny1,Nz1,tau1){}

double compareArrayHostToDevice(double *h, double *d,int size,char *legend)
{
	double h_d[8*CellExtent*CellExtent*CellExtent],t;


	MemoryCopy(h_d,d,size,DEVICE_TO_HOST);

	t = compare(h,h_d,size/sizeof(double),legend,TOLERANCE);

	return t;
}

GPUCell<Particle,dims>* copyCellToDevice( vector<GPUCell<Particle,dims> > * all_cells)
{
	GPUCell<Particle,dims> *h_src,*d_dst;//,*h_ctrl;
	int err1,err2,err3,err4,err5,err6,err7,err8,err9,err10;
	int err11,err12,err13,err14,err15,err16,err17,err18,err19,err20;
	int err21,err22,err23,err24,err25;


	h_src = new GPUCell<Particle,dims>;
	//h_ctrl = new GPUCell<Particle>;

	h_src->number_of_particles = Cell<Particle,dims>::number_of_particles;
	h_src->mesh.x = Cell<Particle,dims>::mesh.x;
	h_src->mesh.y = Cell<Particle,dims>::mesh.y;

	h_src->mesh.setZ(Cell<Particle,dims>::mesh.z());
	h_src->hstep.x = Cell<Particle,dims>::hstep.x;
	h_src->hstep.y = Cell<Particle,dims>::hstep.y;
	h_src->hstep.setZ(Cell<Particle,dims>::hstep.z());
	h_src->cnum.x  = Cell<Particle,dims>::cnum.x;
	h_src->cnum.setZ(Cell<Particle,dims>::cnum.z());
	h_src->cnum.y = Cell<Particle,dims>::cnum.y;
	h_src->xmin.x = Cell<Particle,dims>::xmin.x;
	h_src->xmin.y = Cell<Particle,dims>::xmin.y;
	h_src->xmin.setZ(Cell<Particle,dims>::xmin.z());
	h_src->xmax.x = Cell<Particle,dims>::xmax.x;
	h_src->xmax.y = Cell<Particle,dims>::xmax.y;
	h_src->xmax.setZ(Cell<Particle,dims>::xmax.z());
	h_src->tau = Cell<Particle,dims>::tau;
	h_src->jmp = Cell<Particle,dims>::jmp;
	h_src->d_ctrlParticles = Cell<Particle,dims>::d_ctrlParticles;

	h_src->busyParticleArray = Cell<Particle,dims>::busyParticleArray;

	h_src->beam_boundary = Cell<Particle,dims>::beam_boundary;
	h_src->tgt = Cell<Particle,dims>::tgt;

	//cudaPrintfInit();
	MemoryAllocate((void**)&(h_src->doubParticleArray),sizeof(Particle)*MAX_particles_per_cell);
	err1 = getLastError();



	MemorySet(h_src->doubParticleArray,0,sizeof(Particle)*MAX_particles_per_cell);
	err2 = getLastError();

	//testKernelBefore<<<1,1>>>(h_src->doubParticleArray,50,1);
	//cudaThreadSynchronize();


	MemoryCopy(h_src->doubParticleArray,Cell<Particle,dims>::doubParticleArray,
			   sizeof(Particle)*MAX_particles_per_cell,HOST_TO_DEVICE);
	err3 = getLastError();

	//		(double *)h_src->doubParticleArray,sizeof(Particle)*MAX_particles_per_cell,"part");
	//printf("after copy %e\n",this->ParticleArrayRead(50,1));
	//dbgPrintGPUParticleAttribute(d_dst,50,1," IN_COPY0 " );
	//testKernelBefore<<<1,1>>>(h_src->doubParticleArray,50,1);
	//cudaPrintfDisplay(stdout, true);
	//cudaPrintfEnd();

	MemoryAllocate((void**)&(h_src->Jx),sizeof(CellDouble));
	err4 = getLastError();

	MemoryCopy(h_src->Jx,Cell<Particle,dims>::Jx,sizeof(CellDouble),HOST_TO_DEVICE);
	err5 = getLastError();

	//compareArrayHostToDevice((double *)Cell<Particle,dims>::Jx,(double *)h_src->Jx,sizeof(CellDouble),"Jx");

	MemoryAllocate((void**)&(h_src->Jy),sizeof(CellDouble));
	err6 = getLastError();

	MemoryCopy(h_src->Jy,Cell<Particle,dims>::Jy,sizeof(CellDouble),HOST_TO_DEVICE);
	//compareArrayHostToDevice((double *)Cell<Particle,dims>::Jy,(double *)h_src->Jy,sizeof(CellDouble),"Jy");
	err7 = getLastError();


	MemoryAllocate((void**)&(h_src->Jz),sizeof(CellDouble));
	err8 = getLastError();

	MemoryCopy(h_src->Jz,Cell<Particle,dims>::Jz,sizeof(CellDouble),HOST_TO_DEVICE);
	err9 = getLastError();

	//compareArrayHostToDevice((double *)Cell<Particle,dims>::Jz,(double *)h_src->Jz,sizeof(CellDouble),"Jz");

	MemoryAllocate((void**)&(h_src->Ex),sizeof(CellDouble));
	err10 = getLastError();

	MemoryCopy(h_src->Ex,Cell<Particle,dims>::Ex,sizeof(CellDouble),HOST_TO_DEVICE);
	err11 = getLastError();

	//compareArrayHostToDevice((double *)Cell<Particle,dims>::Ex,(double *)h_src->Ex,sizeof(CellDouble),"Ex");

	MemoryAllocate((void**)&(h_src->Ey),sizeof(CellDouble));
	err12 = getLastError();

	MemoryCopy(h_src->Ey,Cell<Particle,dims>::Ey,sizeof(CellDouble),HOST_TO_DEVICE);
	err13 = getLastError();

	//compareArrayHostToDevice((double *)Cell<Particle,dims>::Ey,(double *)h_src->Ey,sizeof(CellDouble),"Ey");

	MemoryAllocate((void**)&(h_src->Ez),sizeof(CellDouble));
	err14 = getLastError();

	MemoryCopy(h_src->Ez,Cell<Particle,dims>::Ez,sizeof(CellDouble),HOST_TO_DEVICE);
	err15 = getLastError();

	//compareArrayHostToDevice((double *)Cell<Particle,dims>::Ez,(double *)h_src->Ez,sizeof(CellDouble),"Ez");

	MemoryAllocate((void**)&(h_src->Hx),sizeof(CellDouble));
	err16 = getLastError();

	MemoryCopy(h_src->Hx,Cell<Particle,dims>::Hx,sizeof(CellDouble),HOST_TO_DEVICE);
	err17 = getLastError();

	//compareArrayHostToDevice((double *)Cell<Particle,dims>::Hx,(double *)h_src->Hx,sizeof(CellDouble),"Hx");

	MemoryAllocate((void**)&(h_src->Hy),sizeof(CellDouble));
	err18 = getLastError();

	MemoryCopy(h_src->Hy,Cell<Particle,dims>::Hy,sizeof(CellDouble),HOST_TO_DEVICE);
	err19 = getLastError();

	//compareArrayHostToDevice((double *)Cell<Particle,dims>::Hy,(double *)h_src->Hy,sizeof(CellDouble),"Hy");

	MemoryAllocate((void**)&(h_src->Hz),sizeof(CellDouble));
	err20 = getLastError();

	MemoryCopy(h_src->Hz,Cell<Particle,dims>::Hz,sizeof(CellDouble),HOST_TO_DEVICE);
	err21 = getLastError();

	//compareArrayHostToDevice((double *)Cell<Particle,dims>::Hz,(double *)h_src->Hz,sizeof(CellDouble),"Hz");

	MemoryAllocate((void**)&(h_src->Rho),sizeof(CellDouble));
	err22 = getLastError();

	MemoryCopy(h_src->Rho,Cell<Particle,dims>::Rho,sizeof(CellDouble),HOST_TO_DEVICE);
	err23 = getLastError();

	//compareArrayHostToDevice((double *)Cell<Particle,dims>::Rho,(double *)h_src->Rho,sizeof(CellDouble),"Rho");

	//memcpy((unsigned char *)dst.Jx,(unsigned char *)src.Jx,sizeof(CellDouble));
	//printf("i %d l %d k %d q_m %15.5e \n",h_src->i,h_src->k,h_src->l,Cell<Particle,dims>::ParticleArrayRead(0,7));

    MemoryAllocate((void**)&d_dst,sizeof(GPUCell<Particle,dims>));
	err24 = getLastError();



    MemoryCopy(d_dst,h_src,sizeof(GPUCell<Particle,dims>),HOST_TO_DEVICE);
	err25 = getLastError();

	if(
			(err1 != 0) ||
			(err2 != 0) ||
			(err3 != 0) ||
			(err4 != 0) ||
			(err5 != 0) ||
			(err6 != 0) ||
			(err7 != 0) ||
			(err8 != 0) ||
			(err9 != 0) ||
			(err10 != 0) ||
			(err11 != 0) ||
			(err12 != 0) ||
			(err13 != 0) ||
			(err14 != 0) ||
			(err15 != 0) ||
			(err16 != 0) ||
			(err17 != 0) ||
			(err18 != 0) ||
			(err19 != 0) ||
			(err20 != 0) ||
			(err21 != 0) ||
			(err22 != 0) ||
			(err23 != 0) ||
			(err24 != 0) ||
			(err25 != 0)
	  )
	{
		//int qq = 0;
	}

 //   cudaMemcpy(h_ctrl,d_dst,sizeof(Cell<Particle,dims>),cudaMemcpyDeviceToHost);
  //      dbgPrintGPUParticleAttribute(d_dst,50,1," IN_COPY " );
  //  cudaMemcpy(Cell<Particle,dims>::doubParticleArray,h_ctrl->doubParticleArray,
    //			   sizeof(Particle)*MAX_particles_per_cell,cudaMemcpyDeviceToHost);

   // printf("before copy return  %e\n",this->ParticleArrayRead(50,1));
   // dbgPrintGPUParticleAttribute(d_dst,50,1," IN_COPY " );


     return d_dst;
}

void copyCellFromDevice(GPUCell<Particle,dims>* d_src,GPUCell<Particle,dims>* h_dst,string where,int nt)
{
	static GPUCell<Particle,dims> *h_copy_of_d_src;
	static int first = 1;
	int code;


	if(first == 1)
	{
	   first = 0;
	   h_copy_of_d_src = new GPUCell<Particle,dims>;
	   h_copy_of_d_src->Init();

	}

	int err = getLastError();
	if(err != 0)
		{
			 printf(" copyCellFromDevice enter %d %s \n ",err,getErrorString(err));
			 exit(0);
		}


	ThreadSynchronize();

	err = MemoryCopy(h_copy_of_d_src,d_src,sizeof(GPUCell<Particle,dims>),DEVICE_TO_HOST);
	if(err != 0)
	{
		 printf(" copyCellFromDevice1 %d %s \n ",err,getErrorString(err));
		 exit(0);
	}
	//printf("Cuda error: %d: %s.\n", code,getErrorString((cudaError_t) code));
    if(h_copy_of_d_src->number_of_particles < 0 || h_copy_of_d_src->number_of_particles > MAX_particles_per_cell)
    {
    	int qq = 0;
    }
	//code = MemoryCopy(h_dst,&h_copy_of_d_src,sizeof(GPUCell<Particle,dims>),cudaMemcpyHostToHost);
#ifdef COPY_CELLS_MEMORY_PRINTS
	printf("step %d %s number of particles %5d %3d %3d %d \n",nt,where,h_copy_of_d_src->i,h_copy_of_d_src->l,h_copy_of_d_src->k,h_copy_of_d_src->number_of_particles);
#endif



	h_dst->number_of_particles = h_copy_of_d_src->number_of_particles;

	code = MemoryCopy(h_dst->doubParticleArray,h_copy_of_d_src->doubParticleArray,
			   sizeof(Particle)*MAX_particles_per_cell,DEVICE_TO_HOST);
	if(code != 0)
	{
		 printf(" copyCellFromDevice3 %d %s \n ",code,getErrorString(code));
		 exit(0);
	}

	//printf("i %d l %d k %d q_m %15.5e \n",h_dst->i,h_dst->k,h_dst->l,h_dst->ParticleArrayRead(50,1));


	//h_dst->Jx = new CellDouble;
	//compareArrayHostToDevice((double *)h_dst->Jx,(double *)(h_copy_of_d_src.Jx),sizeof(CellDouble),"TEST");

	code = MemoryCopy(h_dst->Jx,h_copy_of_d_src->Jx,sizeof(CellDouble),DEVICE_TO_HOST);
	if(code != 0)
	{
		 printf(" copyCellFromDevice4 %d \n ",code);
		 exit(0);
	}



	//h_dst->Jy = new CellDouble;
	code = MemoryCopy(h_dst->Jy,h_copy_of_d_src->Jy,sizeof(CellDouble),DEVICE_TO_HOST);
	if(code != 0)
	{
		 printf(" copyCellFromDevice5 %d \n ",code);
		 exit(0);
	}

	//h_dst->Jz = new CellDouble;
	code = MemoryCopy(h_dst->Jz,h_copy_of_d_src->Jz,sizeof(CellDouble),DEVICE_TO_HOST);
	if(code != 0)
	{
		 printf(" copyCellFromDevice6 %d \n ",code);
		 exit(0);
	}

	//h_dst->Ex = new CellDouble;
	code = MemoryCopy(h_dst->Ex,h_copy_of_d_src->Ex,sizeof(CellDouble),DEVICE_TO_HOST);
	//h_dst->Ey = new CellDouble;
	code = MemoryCopy(h_dst->Ey,h_copy_of_d_src->Ey,sizeof(CellDouble),DEVICE_TO_HOST);
	//h_dst->Ez = new CellDouble;
	code = MemoryCopy(h_dst->Ez,h_copy_of_d_src->Ez,sizeof(CellDouble),DEVICE_TO_HOST);
	//h_dst->Hx = new CellDouble;
	code = MemoryCopy(h_dst->Hx,h_copy_of_d_src->Hx,sizeof(CellDouble),DEVICE_TO_HOST);
	//h_dst->Hy = new CellDouble;
	code = MemoryCopy(h_dst->Hy,h_copy_of_d_src->Hy,sizeof(CellDouble),DEVICE_TO_HOST);
	//h_dst->Hz = new CellDouble;
	code = MemoryCopy(h_dst->Hz,h_copy_of_d_src->Hz,sizeof(CellDouble),DEVICE_TO_HOST);
	//h_dst->Rho = new CellDouble;
	code = MemoryCopy(h_dst->Rho,h_copy_of_d_src->Rho,sizeof(CellDouble),DEVICE_TO_HOST);
	if(code != 0)
	{
		 printf(" copyCellFromDevice10 %d \n ",code);
		 exit(0);
	}


	//memcpy((unsigned char *)dst.Jx,(unsigned char *)src.Jx,sizeof(CellDouble));
	//printf("i %d l %d k %d q_m %15.5e \n",h_dst->i,h_dst->k,h_dst->l,h_dst->ParticleArrayRead(50,1));

    //cudaMalloc(&d_dst,sizeof(Cell<Particle,dims>));


    //cudaMemcpy(h_dst,d_src,sizeof(Cell<Particle,dims>),cudaMemcpyDeviceToHost);
   // cudaMemcpy(h_ctrl,d_dst,sizeof(Cell<Particle,dims>),cudaMemcpyDeviceToHost);

}

GPUCell<Particle,dims>* allocateCopyCellFromDevice()
{
	GPUCell<Particle,dims> *h_dst,*h_copy_of_d_src;
	//static int first = 1;
	int code;


	   h_dst = new GPUCell<Particle,dims>;
	//h_ctrl = new GPUCell<Particle,dims>;

	h_dst->doubParticleArray = (double*)malloc(sizeof(Particle)*MAX_particles_per_cell);

	h_dst->Jx = new CellDouble;
	h_dst->Jy = new CellDouble;
	h_dst->Jz = new CellDouble;

	h_dst->Ex = new CellDouble;
	h_dst->Ey = new CellDouble;
	h_dst->Ez = new CellDouble;
	h_dst->Hx = new CellDouble;
	h_dst->Hy = new CellDouble;
	h_dst->Hz = new CellDouble;
	h_dst->Rho = new CellDouble;

    return h_dst;
}

void freeCopyCellFromDevice(GPUCell<Particle,dims> *h_dst)
{
	//static int first = 1;
	int code;


	//h_ctrl = new GPUCell<Particle,dims>;

	free(h_dst->doubParticleArray);// = (double*)malloc(sizeof(Particle)*MAX_particles_per_cell);

	delete (h_dst->Jx);// = new CellDouble;
	delete (h_dst->Jy);// = new CellDouble;
	delete (h_dst->Jz);// = new CellDouble;

	delete (h_dst->Ex);// = new CellDouble;
	delete (h_dst->Ey);// = new CellDouble;
	delete (h_dst->Ez);// = new CellDouble;
	delete (h_dst->Hx);// = new CellDouble;
	delete (h_dst->Hy);// = new CellDouble;
	delete (h_dst->Hz);// = new CellDouble;
	delete (h_dst->Rho);// = new CellDouble;

	delete h_dst;
}


#ifdef __CUDACC__
__host__
#endif
void printFileCellParticles(FILE *f,GPUCell<Particle,dims> *h_copy_of_d_src)
{
	Particle p;
	int sorts[3] = {0,0,0};


    fprintf(f,"(%3d,%3d,%3d) ========================================================================================== \n",
    		this->cnum.x,
    		this->cnum.y,
    		this->cnum.z()
    		);

	for(int i = 0;i < h_copy_of_d_src->number_of_particles;i++)
	{
		h_copy_of_d_src->readParticleFromSurfaceDevice(i,&p);
		fprintf(f,"(%3d,%3d,%3d) i %3d sort %d FN %10d c  pointInCell %d %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e %10.3e %10.3e\n",
				this->cnum.x,
				this->cnum.y,
				this->cnum.z(),
				i,(int)p.sort,p.fortran_number,this->isPointInCell(p.GetX()),
		        p.X.x,p.X.y,p.X.z(),p.pu ,p.pv,p.pw,p.q_m,p.m);
		sorts[(int)p.sort] += 1;
	}
	fprintf(f,"ions %3d electrons %3d beam %3d \n",sorts[0],sorts[1],sorts[2]);
}

double compareToCell(Cell<Particle,dims> & d_src)
{

	return Cell<Particle,dims>::compareToCell(d_src);
}

void addParticlesToCellOnDevice(GPUCell<Particle,dims>* d_src,GPUCell<Particle,dims>* h_dst,char *where,int nt)
{
	static GPUCell<Particle,dims> *h_copy_of_d_src;
	static int first = 1;
	int code;


	if(first == 1)
	{
	   first = 0;
	   h_copy_of_d_src = new GPUCell<Particle,dims>;
	   h_copy_of_d_src->Init();

	}

	int err = getLastError();
	if(err != 0)
		{
			 printf(" copyCellFromDevice enter %d %s \n ",err,getErrorString(err));
			 exit(0);
		}


	ThreadSynchronize();

	err = MemoryCopy(h_copy_of_d_src,d_src,sizeof(GPUCell<Particle,dims>),DEVICE_TO_HOST);
	if(err != 0)
	{
		 printf(" copyCellFromDevice1 %d %s \n ",err,getErrorString(err));
		 exit(0);
	}
	//printf("Cuda error: %d: %s.\n", code,getErrorString((cudaError_t) code));
    if(h_copy_of_d_src->number_of_particles < 0 || h_copy_of_d_src->number_of_particles > MAX_particles_per_cell)
    {
    	int qq = 0;
    }
	//code = MemoryCopy(h_dst,&h_copy_of_d_src,sizeof(GPUCell<Particle,dims>),cudaMemcpyHostToHost);
#ifdef COPY_CELLS_MEMORY_PRINTS
	printf("step %d %s number of particles %5d %3d %3d %d \n",nt,where,h_copy_of_d_src->i,h_copy_of_d_src->l,h_copy_of_d_src->k,h_copy_of_d_src->number_of_particles);
#endif



	h_dst->number_of_particles = h_copy_of_d_src->number_of_particles;

	code = MemoryCopy(h_dst->doubParticleArray,h_copy_of_d_src->doubParticleArray,
			   sizeof(Particle)*MAX_particles_per_cell,DEVICE_TO_HOST);
	if(code != 0)
	{
		 printf(" copyCellFromDevice3 %d %s \n ",code,getErrorString(code));
		 exit(0);
	}

	//printf("i %d l %d k %d q_m %15.5e \n",h_dst->i,h_dst->k,h_dst->l,h_dst->ParticleArrayRead(50,1));


	//h_dst->Jx = new CellDouble;
	//compareArrayHostToDevice((double *)h_dst->Jx,(double *)(h_copy_of_d_src.Jx),sizeof(CellDouble),"TEST");

	code = MemoryCopy(h_dst->Jx,h_copy_of_d_src->Jx,sizeof(CellDouble),DEVICE_TO_HOST);
	if(code != 0)
	{
		 printf(" copyCellFromDevice4 %d \n ",code);
		 exit(0);
	}



	//h_dst->Jy = new CellDouble;
	code = MemoryCopy(h_dst->Jy,h_copy_of_d_src->Jy,sizeof(CellDouble),DEVICE_TO_HOST);
	if(code != 0)
	{
		 printf(" copyCellFromDevice5 %d \n ",code);
		 exit(0);
	}

	//h_dst->Jz = new CellDouble;
	code = MemoryCopy(h_dst->Jz,h_copy_of_d_src->Jz,sizeof(CellDouble),DEVICE_TO_HOST);
	if(code != 0)
	{
		 printf(" copyCellFromDevice6 %d \n ",code);
		 exit(0);
	}

	//h_dst->Ex = new CellDouble;
	code = MemoryCopy(h_dst->Ex,h_copy_of_d_src->Ex,sizeof(CellDouble),DEVICE_TO_HOST);
	//h_dst->Ey = new CellDouble;
	code = MemoryCopy(h_dst->Ey,h_copy_of_d_src->Ey,sizeof(CellDouble),DEVICE_TO_HOST);
	//h_dst->Ez = new CellDouble;
	code = MemoryCopy(h_dst->Ez,h_copy_of_d_src->Ez,sizeof(CellDouble),DEVICE_TO_HOST);
	//h_dst->Hx = new CellDouble;
	code = MemoryCopy(h_dst->Hx,h_copy_of_d_src->Hx,sizeof(CellDouble),DEVICE_TO_HOST);
	//h_dst->Hy = new CellDouble;
	code = MemoryCopy(h_dst->Hy,h_copy_of_d_src->Hy,sizeof(CellDouble),DEVICE_TO_HOST);
	//h_dst->Hz = new CellDouble;
	code = MemoryCopy(h_dst->Hz,h_copy_of_d_src->Hz,sizeof(CellDouble),DEVICE_TO_HOST);
	//h_dst->Rho = new CellDouble;
	code = MemoryCopy(h_dst->Rho,h_copy_of_d_src->Rho,sizeof(CellDouble),DEVICE_TO_HOST);
	if(code != 0)
	{
		 printf(" copyCellFromDevice10 %d \n ",code);
		 exit(0);
	}


	//memcpy((unsigned char *)dst.Jx,(unsigned char *)src.Jx,sizeof(CellDouble));
	//printf("i %d l %d k %d q_m %15.5e \n",h_dst->i,h_dst->k,h_dst->l,h_dst->ParticleArrayRead(50,1));

    //cudaMalloc(&d_dst,sizeof(Cell<Particle,dims>));


    //cudaMemcpy(h_dst,d_src,sizeof(Cell<Particle,dims>),cudaMemcpyDeviceToHost);
   // cudaMemcpy(h_ctrl,d_dst,sizeof(Cell<Particle,dims>),cudaMemcpyDeviceToHost);

}





};




#endif /* GPUCELL_H_ */
