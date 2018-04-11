#include "xeon_plasma.h"
#include <stdlib.h>
#include "mpi_shortcut.h"

int main(int argc,char*argv[])
{
      XeonPlasma<GPUCell,DIMENSIONS> *plasma;

      InitMPI(argc,argv);

#ifdef __CUDACC__
      size_t sizeP;

      printf("oarticle size %d %d \n",sizeof(Particle),sizeof(Particle)/sizeof(double));
      cudaDeviceGetLimit(&sizeP,cudaLimitPrintfFifoSize);

      printf("printf default limit %d \n",sizeP/1024/1024);

      sizeP *= 10;
      sizeP *= 10;
      sizeP *= 10;
      sizeP *= 10;
      cudaDeviceSetLimit(cudaLimitPrintfFifoSize, sizeP);

      cudaDeviceGetLimit(&sizeP,cudaLimitPrintfFifoSize);



      printf("printf limit set to %d \n",sizeP/1024/1024);
#endif

      int beam_plasma,start_from_file;
      double tex0,tey0,tez0,Tb,rimp,rbd,ni,lx,ly,lz,tau,B0,bx,by,bz,pl_y,pl_z;
      int lp,nx,ny,nz,np,total_steps,minor_steps,phase;

   readParameterFile(&beam_plasma,
		             &start_from_file,
		             &tex0,&tey0,&tez0,
                     &Tb,&rimp,&rbd,&ni,&lx,&ly,&lz,&lp,&nx,&ny,&nz,
	                 &tau,&B0,&np,&bx,&by,&bz,&pl_y,&pl_z,
	                 &total_steps,&minor_steps,&phase);



   int err = SetDevice(0);
   
      printf("err %d \n",err);
//   int nx0 = 100,ny0 = 4;
//   double hx = 1.1424/nx0,hy = 0.2/ny;
//
//   lx = nx*hx;
//   ly = ny*hy;
   //plasma = new GPUPlasma<GPUCell>(100,4,4,1.2566,0.05,0.05,1.0,100,1.0,0.001);
   plasma = new XeonPlasma<GPUCell,DIMENSIONS>(nx,ny,nz,lx,ly,lz,np,1.0,tau,
		    (DIMENSIONS == 3),B0,80000,
		    tex0,tey0,tez0,Tb,rimp,rbd,ni,bx,by,bz,
		    beam_plasma,start_from_file
		    );

   plasma->SetNumberOfSteps(total_steps,minor_steps,phase);
   plasma->SetPlasmaSize(pl_y,pl_z);

   plasma->Initialize();




   double t = plasma->compareCPUtoGPU();
   printf("----------------------------------------------------------- plasma check before move %.5f\n",t);
   size_t m_free,m_total;

   GetDeviceMemory(&m_free,&m_total);
   struct sysinfo info;
   int start_nt;

   if(START_STEP_NUMBER > 0)
   {
	   start_nt = START_STEP_NUMBER;
   }
   else
   {
	   start_nt = 1;
   }
//   plasma->Diagnose(START_STEP_NUMBER-1);

//   plasma->Bx0 = 0.1;

   for(int nt = START_STEP_NUMBER;nt <= total_steps;nt++)
   {
	   GetDeviceMemory(&m_free,&m_total);
	   sysinfo(&info);
#ifdef MEMORY_PRINTS
       printf("before Step  %10d CPU memory free %10u GPU memory total %10d free %10d\n",
    		   nt,info.freeram/1024/1024,m_total/1024/1024,m_free/1024/1024);
#endif

       plasma->Step(nt);
//       exit(0);
//       plasma->BeamInput(nt);

//       plasma->Diagnose(nt);


       GetDeviceMemory(&m_free,&m_total);
       sysinfo(&info);
#ifdef MEMORY_PRINTS
       printf("after  Step  %10d CPU memory free %10u GPU memory total %10d free %10d\n",
    		   nt,info.freeram/1024/1024,m_total/1024/1024,m_free/1024/1024);
#endif
   }
   exit(0);
   t = plasma->compareCPUtoGPU();
   printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ plasma check after move %.5f\n",t);

   delete plasma;

   CloseMPI();

   return 0;
}
