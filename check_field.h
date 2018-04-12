
#ifndef CHECK_FIELD_H
#define CHECK_FIELD_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ON_DEVICE_CHECK 1
#define CHECK_FIELD_TOLERANCE 1e-4


double checkPointArray(double* a, double* dbg_a,double tol,
		               int size,int silent,int dbg_device_flag,
		               char *name,int num,int nt)
{
	double t = 0.0,*dub_a;
	FILE *f;
	char str[100];
	int wrong = 0;



	sprintf(str,"fch_%s_num%03d_nt%05d.log",name,num,nt);

	if((f = fopen(str,"wt")) == NULL) return 1.0;

	if(dbg_device_flag == 1)
	{
		dub_a = (double *)malloc(size*sizeof(double));

		MemoryCopy(dub_a,dbg_a,size*sizeof(double),DEVICE_TO_HOST);
	}
	else{
		dub_a = dbg_a;
	}

	for(int n = 0;n < size;n++)
	{
		double ta = a[n];
		double t_dbg_a = dub_a[n],q;
		q = (ta - t_dbg_a);
		t =  q >= 0 ? q : -q ;

		if(t >= tol)
		{
           if(silent == 0)
           {
        	   fprintf(f,"%10d correct %15.5e dubious %15.5e diff %15.5e \n",
        		        n,          a[n],       dub_a[n],      t
        		   );
           }
           wrong++;
		}
	}
	fclose(f);

	return (1.0 - ((double)wrong)/size);
}


double printArray	(int num,int nt,char *name,double* a)
	{
	    Cell<Particle,DIMENSIONS> c = (*AllCells)[0];
	    int wrong = 0;
	    double diff = 0.0;
	    char fname[100];
	    FILE *f;


        sprintf(fname,"%s_at%05d_nt_%08d.dat",name,num,nt);
        if((f = fopen(fname,"wt")) == NULL) return 0.0;

	    for(int n = 0;n < mesh.size2();n++)
	    {
		       int3 i = c.getCellTripletNumber(n);
		       fprintf(f,"n %5d i %3d l %3d k %3d %15.5e \n",
				   n,i.x+1,i.y+1,i.z+1,a[n]);
	    }
	    fclose(f);

	    return 1.0;
	}

//double CheckArray	(double* a, double* dbg_a,FILE *f)
//	{
//	    Cell<Particle,DIMENSIONS> c = (*AllCells)[0];
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



double CheckArray	(double* a, double* dbg_a)
	{
	    Cell<Particle,DIMENSIONS> c = (*AllCells)[0];
	    int wrong = 0;
	    double diff = 0.0;
#ifdef CHECK_ARRAY_DETAIL_PRINTS
	    puts("begin array checking2=============================");
#endif
	    for(int n = 0;n < mesh.size2();n++)
	    {
//	        double t  = a[n];
//		    double dt = dbg_a[n];
            diff += pow(a[n] - dbg_a[n],2.0);

	        if(fabs(a[n] - dbg_a[n]) > TOLERANCE)
		    {

		       int3 i = c.getCellTripletNumber(n);
#ifdef CHECK_ARRAY_DETAIL_PRINTS
		       printf("n %5d i %3d l %3d k %3d %15.5e dbg %15.5e diff %15.5e wrong %10d \n",
				   n,i.x+1,i.y+1,i.z+1,a[n],dbg_a[n],fabs(a[n] - dbg_a[n]),wrong++);
#endif
     		}
	    }
#ifdef CHECK_ARRAY_DETAIL_PRINTS
	    printf("  end array checking============================= %.4f less than %15.5e diff %15.5e \n",
	    		(1.0-((double)wrong/(mesh.size2()))),TOLERANCE,
	    		pow(diff/(mesh.size2()),0.5)
	    	  );
#endif

	    return (1.0-((double)wrong/(mesh.size2())));
	}


double CheckArraySilent	(double* a, double* dbg_a)
	{
	    Cell<Particle,DIMENSIONS> c = (*AllCells)[0];
//	    int wrong = 0;
	    double diff = 0.0;
	   // puts("begin array checking=============================");
	    for(int n = 0;n < (mesh.x + 2)*(mesh.y + 2)*mesh.dimz2();n++)
	    {
//	        double t  = a[n];
//		    double dt = dbg_a[n];
            diff += pow(a[n] - dbg_a[n],2.0);

	        if(fabs(a[n] - dbg_a[n]) > TOLERANCE)
		    {

		       int3 i = c.getCellTripletNumber(n);

//		       printf("n %5d i %3d l %3d k %3d %15.5e dbg %15.5e diff %15.5e wrong %10d \n",
//				   n,i.x+1,i.y+1,i.z+1,a[n],dbg_a[n],fabs(a[n] - dbg_a[n]),wrong++);
     		}
	    }
//	    printf("  end array checking============================= %.2f less than %15.5e diff %15.5e \n",
//	    		(1.0-((double)wrong/(mesh.size2()))),TOLERANCE,
//	    		pow(diff/(mesh.size2()),0.5)
//	    	  );

	    return pow(diff/mesh.size2(),0.5);
	}

double CheckGPUArraySilent	(double* a, double* d_a)
	{
	    static double *t;
	    static int f = 1;
	    int err;

	    if(f == 1)
	    {
	    	 t = (double *)malloc(sizeof(double)*mesh.size2());
	    	 f = 0;
	    }

	    MemoryCopy(t,d_a,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
//	    cudaMemcpy(t,d_a,sizeof(double)*mesh.size2(),cudaMemcpyDeviceToHost);
	    err = getLastError();
	    if(err != 0)
	            {
	             	printf("CheckArraySilent err %d %s \n",err,getErrorString(err));
	            	exit(0);
	            }


	   return CheckArraySilent(a,t);
	}

double printGPUArray	(double* d_a,int num,int nt, char *name)
	{
	    static double *t;
	    static int f = 1;
	    int err;

	    if(f == 1)
	    {
	    	 t = (double *)malloc(sizeof(double)*mesh.size2());
	    	 f = 0;
	    }

	    MemoryCopy(t,d_a,sizeof(double)*mesh.size2(),DEVICE_TO_HOST);
//	    cudaMemcpy(t,d_a,sizeof(double)*mesh.size2(),cudaMemcpyDeviceToHost);
	    err = getLastError();
	    if(err != 0)
	            {
	             	printf("printGPUArray err %d %s \n",err,getErrorString(err));
	            	exit(0);
	            }


	   return printArray(num,nt,name,t);
	}



double saveVectorField(double *x,double *y,double *z,
		                int device_flag,int size,char *name,int num,int nt)
{
	FILE *f;
	char str[100];

	sprintf(str,"fch_%s_num%03d_nt%05d.dat",name,num,nt);

	if((f = fopen(str,"wb")) == NULL) return 1.0;

	if(device_flag != ON_DEVICE_CHECK)
	{
	   fwrite(x,sizeof(double),size,f);
	   fwrite(y,sizeof(double),size,f);
	   fwrite(z,sizeof(double),size,f);
	}
	else
	{
		static double *h_x,*h_y,*h_z;
		static int first = 1;

		if(first)
		{
		   h_x = (double *)malloc(size*sizeof(double));
		   h_y = (double *)malloc(size*sizeof(double));
		   h_z = (double *)malloc(size*sizeof(double));
		   first = 0;
		}

		MemoryCopy(h_x,x,size*sizeof(double),DEVICE_TO_HOST);
		MemoryCopy(h_y,y,size*sizeof(double),DEVICE_TO_HOST);
		MemoryCopy(h_z,z,size*sizeof(double),DEVICE_TO_HOST);

		fwrite(h_x,sizeof(double),size,f);
		fwrite(h_y,sizeof(double),size,f);
		fwrite(h_z,sizeof(double),size,f);
	}
	fclose(f);
}

double checkVectorField(double *x,double *y,double *z,
		                int device_flag,int size,string name,int num,int nt)
{
	FILE *f;
	char str[100];

	sprintf(str,"fch_%s_num%03d_nt%05d.dat",name.c_str(),num,nt);

	if((f = fopen(str,"rb")) == NULL) return 1.0;

	static double *h_x,*h_y,*h_z,t_x,t_y,t_z;
	static int first = 1;

	if(first)
	{
	   h_x = (double *)malloc(size*sizeof(double));
	   h_y = (double *)malloc(size*sizeof(double));
	   h_z = (double *)malloc(size*sizeof(double));
	   first = 0;
	}

	fread(h_x,sizeof(double),size,f);
	fread(h_y,sizeof(double),size,f);
	fread(h_z,sizeof(double),size,f);
	
	char cmp[100];
	strcpy(cmp,name.c_str());

	strcat(cmp,"x");
	t_x = checkPointArray(h_x,x,CHECK_FIELD_TOLERANCE,
            size,0,device_flag,
            cmp,num,nt);

	strcpy(cmp,name.c_str());
	strcat(cmp,"y");
	t_y = checkPointArray(h_y,y,CHECK_FIELD_TOLERANCE,
            size,0,device_flag,
            cmp,num,nt);
	strcpy(cmp,name.c_str());
	strcat(cmp,"z");
	t_z = checkPointArray(h_z,z,CHECK_FIELD_TOLERANCE,
            size,0,device_flag,
            cmp,num,nt);

    fclose(f);
    return (t_x+t_y+t_z)/3.0;
}

double saveFields(double *Ex,double *Ey,double *Ez,
		          double *Hx,double *Hy,double *Hz,
		          double *Jx,double *Jy,double *Jz,
		          double *Qx,double *Qy,double *Qz,
                  int device_flag,int size,int num,int nt)
{
//	double t_e,t_h,t_j,t_q;

	saveVectorField(Ex,Ey,Ez,device_flag,size,"EE",num,nt);
	saveVectorField(Hx,Hy,Hz,device_flag,size,"HH",num,nt);
	saveVectorField(Jx,Jy,Jz,device_flag,size,"JJ",num,nt);
	saveVectorField(Qx,Qy,Qz,device_flag,size,"QQ",num,nt);
}

double checkFields(double *Ex,double *Ey,double *Ez,
		          double *Hx,double *Hy,double *Hz,
		          double *Jx,double *Jy,double *Jz,
		          double *Qx,double *Qy,double *Qz,
                  int device_flag,int size,int num,int nt)
{
	double t_e,t_h,t_j,t_q;

#ifndef CHECK_VECTOR_FIELDS
	return 1.0;
#endif

	t_e = checkVectorField(Ex,Ey,Ez,device_flag,size,"EE",num,nt);
	t_h = checkVectorField(Hx,Hy,Hz,device_flag,size,"HH",num,nt);
	t_j = checkVectorField(Jx,Jy,Jz,device_flag,size,"JJ",num,nt);
	t_q = checkVectorField(Qx,Qy,Qz,device_flag,size,"QQ",num,nt);

	printf("CHECK FIELDS point %5d step %5d E %15.4e H %15.4e J %15.4e Q %15.4e \n",num,nt,t_e,t_h,t_j,t_q);

	return (t_e+t_h+t_j+t_q)/4.0;
}



#endif
