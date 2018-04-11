#ifndef PARTICLE_H
#define PARTICLE_H

#include "types.h"

#include <stdio.h>
#include <string>
#include <string.h>
#include <iostream>
#include <math.h>
#include "run_control.h"
#include "control.h"
#include "particle_types.h"
#include "point.h"

#include "plasma/physics.h"

#define GPU_PARTICLE
#define TOLERANCE 1e-15
#define SIZE_TOLERANCE 1e-10

typedef struct CurrentTensorComponent {
	int i11, i12, i13,
	 i21, i22, i23,
	 i31, i32, i33,
	 i41, i42, i43;
	double t[4];
} CurrentTensorComponent;

typedef struct CurrentTensor {
	CurrentTensorComponent Jx,Jy,Jz;
} CurrentTensor;


//typedef struct Tensor CurrentTensor;

typedef char gpu_string[200];


double compare_prints(double *a,double *b,int num,char *legend,double tol,int print_flag)
{
     double t = 0.0;

     for(int i = 0; i < num ;i++)
     {
         if(fabs(a[i] - b[i]) < tol)
         {
            t += 1.0;
#ifdef COMPARE_PRINTS
            if(print_flag == 1) printf(" i %5d a %e b %e diff %e\n",i,a[i],b[i],fabs(a[i] - b[i]));
#endif

         }
         else
         {
#ifdef COMPARE_PRINTS
        	 if(print_flag == 1) printf("WRONG i %5d a %e b %e diff %e\n",i,a[i],b[i],fabs(a[i] - b[i]));
#endif
         }
     }

     if(num > 0) t /= num;
#ifdef COMPARE_PRINTS
     if(print_flag == 1) printf("%30s %.5f\n",legend,t);
#endif
     return t;
}

double compare(double *a,double *b,int num,char *legend,double tol)
{
     return compare_prints(a,b,num,legend,tol,0);
}

int comd(double a,double b)
{
	return (fabs(a - b) < TOLERANCE);
}

class Particle
{
public:  
  
   Point<double,DIMENSIONS,0> X;
   double pu,pv,pw;
   double m,q_m;
   particle_sorts sort;

  // void SetControlSystem(int j,double *c){jmp = j;d_ctrlParticles = c;}
   
#ifdef DEBUG_PLASMA
 //  double3 next_x;
   int fortran_number;

#endif   

#ifdef __CUDACC__
__host__ __device__
#endif
Particle(){}

#ifdef __CUDACC__
__host__ __device__ __forceinline__
#endif
Particle(double x1,double y1, double z1,double u1,double v1,double w1,double m1,double q_m1)
{
   X.x = x1;
   X.y = y1;
   X.setZ(z1);
   pu = u1;
   pv = v1;
   pw = w1;
   m  = m1;
   q_m = q_m1;
}

#ifdef __CUDACC__
__host__ __device__
#endif
  ~Particle(){}

#ifdef __CUDACC__
__host__ __device__
#endif
int getValue(void *result,diagnostics_operations code)
{
	double *p;
	switch(code)
	{
		case DENSITY:         p = (double *)result;
		            		 *p = m;
			                 break;

		case KINETIC_ENERGY: p = (double *)result;
			                *p = (pu * pu + pv * pv + pw * pw) + 1.0;
			                 break;

		case VELOCITY:       double ps = pow(((pu * pu + pv * pv + pw * pw) + 1.0),-0.5);
		                     double3 *d;
		                     d->x = ps * pu;
		                     d->y = ps * pv;
		                     d->z = ps * pw;
		                     result = d;
			                 break;
	}

	return 0;
}
  
#ifdef __CUDACC__
__host__ __device__ __forceinline__
#endif
void Move(double3 E,double3 H,double tau,double *p_control,int jmp_control)
{
    double bx,by,bz,tau1,u,v,w,ps,su,sv,sw,s1,s2,s3,s4,s5,s6,s;
	double sx,sy,sz,x1,y1,z1,pu1,pv1,pw1;

#ifdef ATTRIBUTES_CHECK
	p_control[ParticleAttributePosition(jmp_control,fortran_number,sort,1)] = x;
	p_control[ParticleAttributePosition(jmp_control,fortran_number,sort,2)] = y;
	p_control[ParticleAttributePosition(jmp_control,fortran_number,sort,3)] = z;
	p_control[ParticleAttributePosition(jmp_control,fortran_number,sort,4)] = pu;
	p_control[ParticleAttributePosition(jmp_control,fortran_number,sort,5)] = pv;
	p_control[ParticleAttributePosition(jmp_control,fortran_number,sort,6)] = pw;
#endif

	tau1=q_m*tau*0.5;
    
	pu += tau1*E.x;
	pv += tau1*E.y;
	pw += tau1*E.z;

#ifdef __CUDACC__
	ps = tau1 * rsqrt((pu * pu + pv * pv + pw * pw) * 1. + 1.0);
#else
	ps = tau1 * sqrt((pu * pu + pv * pv + pw * pw) * 1. + 1.0);
#endif

	bx = ps * H.x;
	by = ps * H.y;
	bz = ps * H.z;
	su = pu + pv * bz - pw * by;
	sv = pv + pw * bx - pu * bz;
	sw = pw + pu * by - pv * bx;

	s1 = bx * bx;
	s2 = by * by;
	s3 = bz * bz;
	s4 = bx * by;
	s5 = by * bz;
	s6 = bz * bx;
	s = s1 + 1. + s2 + s3;

	sx = tau1*E.x;
	sy = tau1*E.y;
	sz = tau1*E.z;

	pu1 = ((s1 + 1.) * su + (s4 + bz) * sv + (s6 - by) * sw) / s;
	pv1 = ((s4 - bz) * su + (s2 + 1.) * sv + (s5 + bx) * sw) / s;
	pw1 = ((s6 + by) * su + (s5 - bx) * sv + (s3 + 1.) * sw) / s;

	pu = pu1 + sx;
	pv = pv1 + sy;
	pw = pw1 + sz;
	ps = pu * pu + pv * pv + pw * pw;
	ps = pow(((pu * pu + pv * pv + pw * pw) + 1.0),-0.5);

	u = ps * pu;
	v = ps * pv;
	w = ps * pw;
	x1 = X.x   + tau * u;
	y1 = X.y   + tau * v;
	z1 = X.z() + tau * w;
	

	X.x = x1;
	X.y = y1;
	X.setZ(z1);
	
}

#ifdef __CUDACC__
__host__ __device__
#endif
void Collide(double sect){}

#ifdef __CUDACC__
__host__ __device__ __forceinline__
#endif
   double3 GetX(){double3 d3x; d3x.x = X.x; d3x.y = X.y; d3x.z = X.z(); return d3x;}

//#ifdef __CUDACC__
//__host__ __device__ __forceinline__
//#endif
//   double3 GetV(){double3 d3x; d3x.x = x; d3x.y = y; d3x.z = z; return d3x;}

#ifdef __CUDACC__
__host__ __device__ __forceinline__
#endif
   void    SetX(double3 x1){X.x = x1.x;X.y = x1.y;X.setZ(x1.z);}

#ifdef __CUDACC__
__host__ __device__ __forceinline__
#endif
   double  GetMass(){return m;}

#ifdef __CUDACC__
__host__ __device__ __forceinline__
#endif
   double  GetQ2M(){return q_m;}
   

#ifdef DEBUG_PLASMA



#endif


void Print(FILE* f, int num)
{
   }

#ifdef __CUDACC__
__host__ __device__
#endif
Particle & operator=(Particle & src)
{
	X.x = src.X.x;
	X.y = src.X.y;
	double z = (src.X).z();
	X.setZ(z);
	pu = src.pu;
	pv = src.pv;
	pw = src.pw;
	m   = src.m;
	q_m = src.q_m;
	sort = src.sort;


#ifdef DEBUG_PLASMA
	fortran_number = src.fortran_number;
#endif
	return (*this);
}

};

#endif
