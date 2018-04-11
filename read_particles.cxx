#include<stdio.h>
#include<stdlib.h>

#include <sys/resource.h>
#include <stdint.h>

#include <sys/sysinfo.h>
#include <sys/time.h>

#include "particle_types.h"

//#include "read_particles.h"

#include "archAPI.h"

#include <string.h>

int readBinaryParticleArraysOneSort(
    		  FILE *f,
    		  double **dbg_x,
    		  double **dbg_y,
    		  double **dbg_z,
    		  double **dbg_px,
    		  double **dbg_py,
    		  double **dbg_pz,
    		  double *qq_m,
    		  double **mm,
    		  int nt,
    		  int sort
    		  )
      {
		     double q_m,tp,m;
		     int t;
		     //Cell<Particle> c0 = (*AllCells)[0];
		     int total_particles;
		     int err;
		     long int pos = ftell(f);

		     if((err = ferror(f)) != 0)
		     {
		     	 return err ;
		     }

		     fread(&t,sizeof(int),1,f);
		     if((err = ferror(f)) != 0)
		    	 {
		    	 	 return err ;
		    	 }
		     fread(&tp,sizeof(double),1,f);
		     if((err = ferror(f)) != 0)
		    	 {
		    	 	 return err ;
		    	 }

		     total_particles = (int)tp;
		     fread(&q_m,sizeof(double),1,f);
		     if((err = ferror(f)) != 0)
		    	 {
		    	 	 return err ;
		    	 }

		     fread(&m,sizeof(double),1,f);
		     if((err = ferror(f)) != 0)
		    	 {
		    	 	 return err ;
		    	 }

		     fread(&t,sizeof(int),1,f);
		     if((err = ferror(f)) != 0)
		    	 {
		    	 	 return err ;
		    	 }

	         *dbg_x = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_x,total_particles,nt,"x",sort);

	         *dbg_y = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_y,total_particles,nt,"y",sort);

	         *dbg_z = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_z,total_particles,nt,"z",sort);

	         *dbg_px = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_px,total_particles,nt,"px",sort);

	         *dbg_py = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_py,total_particles,nt,"py",sort);

	         *dbg_pz = (double *)malloc(sizeof(double)*total_particles);

	         *mm     = (double *)malloc(sizeof(double)*total_particles);

	         //debugPrintParticleCharacteristicArray(*dbg_pz,total_particles,nt,"pz",sort);

		 	readFortranBinaryArray(f,*dbg_x);
		 	readFortranBinaryArray(f,*dbg_y);
		 	readFortranBinaryArray(f,*dbg_z);
		 	readFortranBinaryArray(f,*dbg_px);
		 	readFortranBinaryArray(f,*dbg_py);
		 	readFortranBinaryArray(f,*dbg_pz);
		 	readFortranBinaryArray(f,*mm);

		 	*qq_m = q_m;
//		 	*mm   = m;

		 	if((err = ferror(f)) != 0)
            {
	   	 	    return err ;
			}

		 	return total_particles;
      }

int AllocateBinaryParticleArraysOneSort(
    		  double **dbg_x,
    		  double **dbg_y,
    		  double **dbg_z,
    		  double **dbg_px,
    		  double **dbg_py,
    		  double **dbg_pz,
    		  double **m,
    		  int total_particles
    		  )
      {
	         *dbg_x = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_x,total_particles,nt,"x",sort);

	         *dbg_y = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_y,total_particles,nt,"y",sort);

	         *dbg_z = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_z,total_particles,nt,"z",sort);

	         *dbg_px = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_px,total_particles,nt,"px",sort);

	         *dbg_py = (double *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_py,total_particles,nt,"py",sort);

	         *dbg_pz = (double *)malloc(sizeof(double)*total_particles);

	         *m      = (double *)malloc(sizeof(double)*total_particles);

	         //debugPrintParticleCharacteristicArray(*dbg_pz,total_particles,nt,"pz",sort);

		 	return 0;
      }

int AllocateBinaryParticleArraysOneSortFloat(
    		  float **dbg_x,
    		  float **dbg_y,
    		  float **dbg_z,
    		  float **dbg_px,
    		  float **dbg_py,
    		  float **dbg_pz,
    		  int total_particles
    		  )
      {
	         *dbg_x = (float*)malloc(sizeof(float)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_x,total_particles,nt,"x",sort);

	         *dbg_y = (float*)malloc(sizeof(float)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_y,total_particles,nt,"y",sort);

	         *dbg_z = (float *)malloc(sizeof(float)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_z,total_particles,nt,"z",sort);

	         *dbg_px = (float *)malloc(sizeof(float)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_px,total_particles,nt,"px",sort);

	         *dbg_py = (float *)malloc(sizeof(float)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_py,total_particles,nt,"py",sort);

	         *dbg_pz = (float *)malloc(sizeof(double)*total_particles);
	         //debugPrintParticleCharacteristicArray(*dbg_pz,total_particles,nt,"pz",sort);

		 	return 0;
      }


	int readFortranBinaryArray(FILE *f, double* d)
	{
//	    Cell<Particle>  c = (*AllCells)[0];
	    int t,err;//,n;


	    fread(&t,sizeof(int),1,f);
	     if((err = ferror(f)) != 0)
	    	 {
	    	 	 return err ;
	    	 }

	    fread(d,1,t,f);
	     if((err = ferror(f)) != 0)
	    	 {
	    	 	 return err ;
	    	 }

	    fread(&t,sizeof(int),1,f);
	     if((err = ferror(f)) != 0)
	    	 {
	    	 	 return err ;
	    	 }



#ifdef READ_DEBUG_PRINTS
	    for(int i = 1; i <= Nx+2;i++)
	    {
	    	for(int l = 1; l <= Ny+2;l++)
	    	{
	    		for(int k = 1;k <= Nz+2;k++)
	    		{
	    			n = c.getFortranCellNumber(i,l,k);
	    			printf("%5d %5d %5d %25.15e \n",i,l,k,d[n]);
	    		}
	    	}
	    }
#endif

	    	    return t;
	}


	  void InitBinaryParticlesArrays(char *fn,int nt,
			  ParticleArrays *ions,
			  ParticleArrays *electrons,
			  ParticleArrays *beam_electrons,
			  int Nx,int Ny,int Nz,
			  int beam_plasma
			  )
	  {
	     FILE *f;
	     double *buf;
	     long int pos;

	     buf = (double *)malloc(sizeof(double)*(Nx+2)*(Ny+2)*(Nz+2));

	     if((f = fopen(fn,"rb")) == NULL) return;
	     struct sysinfo info;

	     pos = ftell(f);

	     sysinfo(&info);
	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
	     readFortranBinaryArray(f,buf); //ex
	     pos = ftell(f);
	     sysinfo(&info);
	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
	     readFortranBinaryArray(f,buf); //ey
	     pos = ftell(f);
	     sysinfo(&info);
	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);

	     readFortranBinaryArray(f,buf);  //ez
	     pos = ftell(f);

	     sysinfo(&info);
	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);

	     readFortranBinaryArray(f,buf);  //hx
	     pos = ftell(f);

	     sysinfo(&info);
	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);

	     readFortranBinaryArray(f,buf);   //hy
	     pos = ftell(f);

	     readFortranBinaryArray(f,buf);   //hz
	     pos = ftell(f);

	     readFortranBinaryArray(f,buf);   //jx
	     pos = ftell(f);

	     readFortranBinaryArray(f,buf);   //jy
	     pos = ftell(f);

	     readFortranBinaryArray(f,buf);   //jz
	     pos = ftell(f);

//	     readFortranBinaryArray(f,buf);
//
//	     readFortranBinaryArray(f,buf);
//	     readFortranBinaryArray(f,buf);
//	     readFortranBinaryArray(f,buf);
	     int err;
	     err = ferror(f);

	     ions->total = readBinaryParticleArraysOneSort(f,&(ions->dbg_x), &(ions->dbg_y) ,&(ions->dbg_z),
	    		                                         &(ions->dbg_px),&(ions->dbg_py),&(ions->dbg_pz),
	    		                                             &(ions->q_m),&(ions->m),nt,ION);

//	     readBinaryParticlesOneSort(f,vp,ION,nt,dbg_ion_x,dbg_ion_y,dbg_ion_z,
//                                               dbg_ion_px,dbg_ion_py,dbg_ion_pz,
//                                               ion_q_m,ion_m,total_ions);

	     sysinfo(&info);
	     printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
	     err = ferror(f);
         if(beam_plasma == 1)
         {
//	     readBinaryParticlesOneSort1(f,vp,PLASMA_ELECTRON,nt);
	        electrons->total = readBinaryParticleArraysOneSort(f,&(electrons->dbg_x),&(electrons->dbg_y), &(electrons->dbg_z),
	    		                                             &(electrons->dbg_px),&(electrons->dbg_py),&(electrons->dbg_pz),
	    		                                             &(electrons->q_m),&(electrons->m),nt,PLASMA_ELECTRON);

            err = ferror(f);
	        sysinfo(&info);
	        printf("before1  %d free %u \n",nt,info.freeram/1024/1024);
            err = ferror(f);
//	     readBinaryParticlesOneSort1(f,vp,BEAM_ELECTRON,nt);
	        beam_electrons->total = readBinaryParticleArraysOneSort(f,&(beam_electrons->dbg_x),&(beam_electrons->dbg_y),&(beam_electrons->dbg_z),
	    		                                             &(beam_electrons->dbg_px),&(beam_electrons->dbg_py),&(beam_electrons->dbg_pz),
	    		                                             &(beam_electrons->q_m),&(beam_electrons->m),nt,BEAM_ELECTRON);

//	     readBinaryParticlesOneSort(f,vp,BEAM_ELECTRON,nt,dbg_beam_x,dbg_beam_y,dbg_beam_z,
//                                               dbg_beam_px,dbg_beam_py,dbg_beam_pz,
//                                               beam_q_m,beam_m,total_beam_electrons);

	        fclose(f);
         }

	     //magf = 1;
	  }

	  void AllocateBinaryParticlesArrays(
			  ParticleArrays *ions,ParticleArrays *electrons,ParticleArrays *beam_electrons)
	  {
	     AllocateBinaryParticleArraysOneSort(&(ions->dbg_x), &(ions->dbg_y) ,&(ions->dbg_z),
	    		                                         &(ions->dbg_px),&(ions->dbg_py),
	    		                                         &(ions->dbg_pz),
	    		                                         &(ions->m),
	    		                                          ions->total);

	     AllocateBinaryParticleArraysOneSort(&(electrons->dbg_x),&(electrons->dbg_y), &(electrons->dbg_z),
	    		                                             &(electrons->dbg_px),&(electrons->dbg_py),&(electrons->dbg_pz),
	    		                                             &(electrons->m),
	    		                                             electrons->total);

	     AllocateBinaryParticleArraysOneSort(&(beam_electrons->dbg_x),&(beam_electrons->dbg_y),&(beam_electrons->dbg_z),
	    		                                             &(beam_electrons->dbg_px),&(beam_electrons->dbg_py),&(beam_electrons->dbg_pz),
	    		                                             &(beam_electrons->m),
	    		                                             beam_electrons->total);

	     //magf = 1;
	  }

	  void AllocateBinaryParticlesArraysFloat(
			  ParticleFloatArrays *ions,ParticleFloatArrays *electrons,ParticleFloatArrays *beam_electrons)
	  {
	     AllocateBinaryParticleArraysOneSortFloat(&(ions->dbg_x), &(ions->dbg_y) ,&(ions->dbg_z),
	    		                                         &(ions->dbg_px),&(ions->dbg_py),&(ions->dbg_pz),
	    		                                          ions->total);


	     AllocateBinaryParticleArraysOneSortFloat(&(electrons->dbg_x),&(electrons->dbg_y), &(electrons->dbg_z),
	    		                                             &(electrons->dbg_px),&(electrons->dbg_py),&(electrons->dbg_pz),
	    		                                             electrons->total);

	     AllocateBinaryParticleArraysOneSortFloat(&(beam_electrons->dbg_x),&(beam_electrons->dbg_y),&(beam_electrons->dbg_z),
	    		                                             &(beam_electrons->dbg_px),&(beam_electrons->dbg_py),&(beam_electrons->dbg_pz),
	    		                                             beam_electrons->total);

	     //magf = 1;
	  }

	  int AllocateDeviceParticleDiagnosticPointers(ParticleFloatArraysGroup **d_pfag,
			                                       ParticleFloatArraysGroup *host_copy_d_pfag,
			                                       ParticleFloatArraysGroup *pfag)
	  {
//		  ParticleFloatArraysGroup host_copy_d_pfag;
		  int i,tot;

		  for(i = 0;i < SORTS;i++)
		  {
			  tot = (*pfag)[i].total;

			  MemoryAllocate((void **)&( (*host_copy_d_pfag)[i].dbg_px),tot*sizeof(float));
			  MemoryAllocate((void **)&( (*host_copy_d_pfag)[i].dbg_py),tot*sizeof(float));
			  MemoryAllocate((void **)&( (*host_copy_d_pfag)[i].dbg_pz),tot*sizeof(float));
			  MemoryAllocate((void **)&( (*host_copy_d_pfag)[i].dbg_x),tot*sizeof(float));
			  MemoryAllocate((void **)&( (*host_copy_d_pfag)[i].dbg_y),tot*sizeof(float));
			  MemoryAllocate((void **)&( (*host_copy_d_pfag)[i].dbg_z),tot*sizeof(float));

			  (*host_copy_d_pfag)[i].total = (*pfag)[i].total;

			  (*host_copy_d_pfag)[i].m = (*pfag)[i].m;
			  (*host_copy_d_pfag)[i].q_m = (*pfag)[i].q_m;
		  }
		  MemoryAllocate((void **)d_pfag,sizeof(ParticleFloatArraysGroup));
		  int err = MemoryCopy(*d_pfag,host_copy_d_pfag,sizeof(ParticleFloatArraysGroup),HOST_TO_DEVICE);

		  return 0;
	  }

	  int CopyParticleDiagnosticPointersFromDevice(ParticleFloatArraysGroup pfag,
			                                       ParticleFloatArraysGroup host_copy_d_pfag)
	  {
//		  ParticleFloatArraysGroup host_copy_d_pfag;
		  int i,tot;

		  for(i = 0;i < SORTS;i++)
		  {
			  tot = pfag[i].total;

			  MemoryCopy(pfag[i].dbg_px,host_copy_d_pfag[i].dbg_px,tot*sizeof(float),DEVICE_TO_HOST);
			  MemoryCopy(pfag[i].dbg_py,host_copy_d_pfag[i].dbg_py,tot*sizeof(float),DEVICE_TO_HOST);
			  MemoryCopy(pfag[i].dbg_pz,host_copy_d_pfag[i].dbg_pz,tot*sizeof(float),DEVICE_TO_HOST);
			  MemoryCopy(pfag[i].dbg_x,host_copy_d_pfag[i].dbg_x,tot*sizeof(float),DEVICE_TO_HOST);
			  MemoryCopy(pfag[i].dbg_y,host_copy_d_pfag[i].dbg_y,tot*sizeof(float),DEVICE_TO_HOST);
			  MemoryCopy(pfag[i].dbg_z,host_copy_d_pfag[i].dbg_z,tot*sizeof(float),DEVICE_TO_HOST);
		  }

		  return 0;
	  }


