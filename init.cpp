#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "rnd.h"

#include "read_particles.h"

#include "run_control.h"



// #define COMPARE_PRINTS

//#define N 160000
//
//double xi[N],yi[N],zi[N],ui[N],vi[N],wi[N];
//double xb[N],yb[N],zb[N],ub[N],vb[N],wb[N];
//double xf[2*N],yf[2*N],zf[2*N],uf[2*N],vf[2*N],wf[2*N];

/* plasma.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"
#include <stdio.h>
// #include "run_control.h"

/* Common Block Declarations */

struct cag05b_1_ {
    doublereal store1, store2;
};
struct cag05b_2_ {
    doublereal normal, gamma;
};

#define cag05b_1 (*(struct cag05b_1_ *) &cag05b_)
#define cag05b_2 (*(struct cag05b_2_ *) &cag05b_)

struct cag05a_1_ {
    integer ix, iy, iz;
};

#define cag05a_1 (*(struct cag05a_1_ *) &cag05a_)

/* Initialized data */

struct {
    doublereal e_1[2];
    } cag05b_ = { 1., -1. };

struct {
    integer e_1[3];
    } cag05a_ = { 1, 255, 25555 };


/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__0 = 0;
static integer c__2 = 2;
static doublereal c_b172 = 1.;
static integer c__12 = 12;
static doublereal c_b417 = 1.1424;
static doublereal c_b419 = .5712;
static doublereal c_b429 = 1.5;
static integer c__20000 = 20000;
static doublereal c_b581 = 0.;
static doublereal c_b582 = .11200000000000002;
static doublereal c_b587 = .14;
static doublereal c_b589 = .8;
static doublereal c_b614 = .001;

/* ------------------------------------------------------ */
/* ������� �.�., */
doublereal g05dde_(doublereal *a, doublereal *b,int dbg_print)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal half = .5;
    static doublereal d__[41] = { 0.,.674489750196082,1.150349380376008,
	    1.534120544352546,1.862731867421652,2.153874694061456,
	    2.417559016236505,2.66006746861746,2.885634912426757,
	    3.097269078198785,3.297193345691964,3.487104104114431,
	    3.668329285121323,3.841930685501911,4.008772594168585,
	    4.169569323349106,4.324919040826046,4.475328424654204,
	    4.621231001499247,4.763001034267814,4.900964207963193,
	    5.035405969463927,5.166578119728753,5.294704084854598,
	    5.419983174916868,5.54259405780294,5.662697617459439,
	    5.780439324478935,5.895951216739571,6.009353565530745,
	    6.120756285971941,6.230260137989044,6.33795775455379,
	    6.443934526538564,6.548269367831731,6.651035379893011,
	    6.752300431407015,6.852127665896068,6.95057594791675,
	    7.047700256664409,7.14355203435219 };

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer n;
    static doublereal t, u, v, w, x;
    extern doublereal wrapg05cae_(doublereal *,int);

    u = cag05b_1.store1;
    for (n = 1; n <= 39; ++n) {
	if (u > half) {
	    goto L40;
	}
	u += u;
/* L20: */
    }
    n = 40;
L40:
    t = d__[n - 1];
    u = wrapg05cae_(&x,dbg_print);
L60:
    w = (d__[n] - t) * u;
    v = w * (w * half + t);
L80:
    u = wrapg05cae_(&x,dbg_print);
    if (v <= u) {
	goto L100;
    }
    v = wrapg05cae_(&x,dbg_print);
    if (u > v) {
	goto L80;
    }
    u = (v - u) / (one - u);
    goto L60;
L100:
    u = (u - v) / (one - v);
    if (u > half) {
	goto L120;
    }
    cag05b_1.store1 = u + u;
    ret_val = *a + *b * (w + t);
    return ret_val;
L120:
    cag05b_1.store1 = u + u - one;
    ret_val = *a - *b * (w + t);
    return ret_val;
} /* g05dde_ */

/* ------------------------------------------------------------------ */
doublereal g05cae_(doublereal *x,int dbg_print)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal ai;
    static integer ii;
    static doublereal ax, ay, az;

    cag05a_1.ix = (cag05a_1.ix - cag05a_1.ix / 177 * 177) * 171 - (
	    cag05a_1.ix / 177 << 1);
    cag05a_1.iy = (cag05a_1.iy - cag05a_1.iy / 176 * 176) * 172 - (
	    cag05a_1.iy / 176 << 1);
    cag05a_1.iz = (cag05a_1.iz - cag05a_1.iz / 178 * 178) * 170 - (
	    cag05a_1.iz / 178 << 1);
    if (cag05a_1.ix < 0) {
	cag05a_1.ix += 30269;
    }
    if (cag05a_1.iy < 0) {
	cag05a_1.iy += 30307;
    }
    if (cag05a_1.iz < 0) {
	cag05a_1.iz += 30323;
    }
    ax = (doublereal) cag05a_1.ix;
    ay = (doublereal) cag05a_1.iy;
    az = (doublereal) cag05a_1.iz;
    ai = ax / 30269. + ay / 30307. + az / 30323.;
    ii = (integer) ai;
    ret_val = ai - ii;
    return ret_val;
} /* g05cae_ */

doublereal wrapg05cae_(doublereal *x,int dbg_print)
{
    static int n = 0;
    double t  = g05cae_(x,dbg_print);
// #ifdef DEBUG_PLASMA
    n++;
    if(dbg_print == 1)
    {
     printf("%10d %25.15e \n",n,t);
    }
// #endif

    return t;
}

doublereal wrapg05dde_(doublereal *a, doublereal *b,int dbg_print)
{
    return g05dde_(a,b,dbg_print);
}

double rnd_uniform(int dbg_print)
{
    doublereal x;

    return (double)wrapg05cae_(&x,dbg_print);
}

double rnd_gaussian(double a,double b,int dbg_print)
{
    return (double)wrapg05dde_(&a,&b,dbg_print);
}


int in_range(double z0,double z,double z1)
{
	return ((z > z0) && (z < z1)) || ((fabs(z - z0) < 1e-13) && (fabs(z - z1) < 1e-13));
}

int InitUniformMaxwellianParticles(int beamf,int jmb,
				   double tex0,double tey0,double tez0,
				   double beam_lx, double beam_ly,double beam_lz,int *jmb_real,
                   double lx,double ly,double lz, int meh,double Tb,double rimp,double rbd,
				   double *xi,double *yi, double *zi,double *ui,double *vi, double *wi,
				   double *xb,double *yb, double *zb,double *ub,double *vb, double *wb,
				   double *xf,double *yf, double *zf,double *uf,double *vf, double *wf
				  )
{
    double x,y,z,vb0,d__1,d__2,d__3,vy,vz,termx,gb0;
    double vf01,vf02,pinv1,pinv2,mfrq = 0.0;
//     double *ux,*uy,*uz;
    double *ux,*uy,*uz;
    double beam_y_max,beam_y_min, beam_sh;
    double beam_z_max,beam_z_min, beam_shz;

    beam_sh = (ly - beam_ly)/2;
    beam_y_max = ly - beam_sh;
    beam_y_min = beam_sh;

    beam_shz   = (lz - beam_lz)/2;
    beam_z_max = lz - beam_shz;
    beam_z_min = beam_shz;

    int j;
    
    ux = (double *)malloc(jmb*sizeof(double));
    uy = (double *)malloc(jmb*sizeof(double));
    uz = (double *)malloc(jmb*sizeof(double));
    
    for (j = 1; j <= jmb; j++) 
    {
	z =   lz * rnd_uniform(0);
	y =   meh * ly + ly * rnd_uniform(0);
	x =   lx * rnd_uniform(0);
	
	xi[j - 1] = x;
	yi[j - 1] = y;
	zi[j - 1] = z;
	ui[j - 1] = 0.0;
	vi[j - 1] = 0.0;
	wi[j - 1] = 0.0;
    }
    
//    parasetrandombeam_();
/*     END RANDOM GENERATOR */
/* ****************** BEAM **************************************** */
    *jmb_real = 0;
	for (j = 1; j <= jmb; j++) 
	{
		double y = yi[j- 1];
		double z = zi[j- 1];
		if((xi[j - 1] < beam_lx) &&
				(y < beam_y_max) && (y > beam_y_min) &&
				in_range(beam_z_min,z,beam_z_max)
		  )
		{
	        xb[*jmb_real] = xi[j - 1];
	        yb[*jmb_real] = yi[j - 1];
	        zb[*jmb_real] = zi[j - 1];
	        vb0       = rnd_gaussian(0.0, Tb*rimp,0);
	        ux[*jmb_real] = vb0 + rimp;
	        uy[*jmb_real] = rnd_gaussian(0.0, Tb*rimp,0);
	        uz[*jmb_real] = rnd_gaussian(0.0, Tb*rimp,0);
#ifdef DEBUG_INITIAL_PARTICLE_PRINTS
            	printf("ion %10d %25.15e %25.15e %25.15e \n",j,xi[j - 1],yi[j - 1],zi[j - 1]);
#endif
	        (*jmb_real)++;
		}
	    
	}
	
	for (j = 1; j <= *jmb_real; j++)
	{
	    double uxt,ubt; 
	    d__1 = ux[j - 1];
	    d__2 = uy[j - 1];
	    d__3 = uz[j - 1];
	    
	    vb0 = sqrt(1.0 - d__1 * d__1 - d__2 * d__2 - d__3 * d__3);
// 	    uxt = ux[j - 1];
// 	    ubt = uxt/vb0;
// 	    ub[j-1] = ubt;
	    ub[j - 1] = ux[j - 1] / vb0;
	    vb[j - 1] = uy[j - 1] / vb0;
	    wb[j - 1] = uz[j - 1] / vb0;
#ifdef DEBUG_INITIAL_PARTICLE_PRINTS
	       printf("beam %10d %25.15e  %25.15e    %25.15e %25.15e %25.15e   %25.15e \n",
	    		         j,  xb[j - 1],yb[j - 1],vb0,    ub[j-1],vb[j - 1],wb[j - 1]);
#endif
	}
/*     MAKING THE RANDOM GENERATOR WORK THE SAME */
//    parasetrandomelectrons_();
/*     END RANDOM GENERATOR */

//     vy = rnd_gaussian(0.0,tey0,1);
//     vz = rnd_gaussian(0.0,tez0,1);
//  
//     termx = rnd_gaussian(0.0,tex0,1); 
	j = 1;
    for (j = 1; j <= jmb;j++) 
    {
    	   if((2*j-1) == 24933)
    	   {
    		   int qq = 0;
    	   }
           xf[2*j-1-1] = xi[j-1];
           yf[2*j-1-1] = yi[j-1];
           zf[2*j-1-1] = zi[j-1];

           xf[2*j-1]   = xi[j-1];
           yf[2*j-1]   = yi[j-1];
           zf[2*j-1]   = zi[j-1];
	   
//         FIRST SETTING TRANVERSE
// razbros v skorostyax
           
           vy=rnd_gaussian(0.0,tey0,0);    
           vz=rnd_gaussian(0.0,tez0,0);        

//          INVERSE CURRENT
           
           termx = rnd_gaussian(0.0,tex0,0);
// 	   printf("termx %15.5e vx %15.5e vy %15.5e \n",termx,vy,vz);
	   
           gb0 = pow(1.0+pow(ub[j-1],2)+pow(vb[j-1],2)+pow(wb[j-1],2),-0.5);
	   
           vb0=ub[j-1]*gb0;
      	   if ((beamf == 1) && ((xi[j - 1] < beam_lx) && (yi[j - 1] < beam_y_max) && (yi[j - 1] > beam_y_min)))
     	   {
        	    vf01=-rbd*vb0+termx;
                vf02=-rbd*vb0-termx;
	       }
       	   else
	       {
      	        vf01=+termx;
      	        vf02=-termx;
	       }
       	   
           pinv1= vf01*pow((1.0-pow(vf01,2)-vy*vy-vz*vz),-0.5); 
           pinv2= vf02*pow((1.0-pow(vf02,2)-vy*vy-vz*vz),-0.5);
            
          vf[2*j-2] =  vy*pow((1.0-pow(vf01,2)-vy*vy-vz*vz),-0.5);
// 	  write(37,127) 2*j-1,vy,vf01,vz,                                                                                                                                                                        
//      +    vy/dsqrt((1d0-vf01**2-vy**2-vz**2))
//           printf("%10d vy %15.5e vf01 %15.5e vz %15.5e vf %15.5e \n",2*j-1,vy,vf01,vz,vy*pow((1.0-pow(vf02,2)-vy*vy-vz*vz),-0.5));
          vf[2*j-1] = -vy*pow((1.0-pow(vf02,2)-vy*vy-vz*vz),-0.5);
	  
          wf[2*j-2] =  vz*pow((1.0-pow(vf01,2)-vy*vy-vz*vz),-0.5);
          wf[2*j-1] = -vz*pow((1.0-pow(vf02,2)-vy*vy-vz*vz),-0.5);

	      uf[2*j-2] = pinv1+0.01*sin(mfrq*2.0*M_PI*xf[2*j-2]/lx);
          uf[2*j-1] = pinv2+0.01*sin(mfrq*2.0*M_PI*xf[2*j-1]/lx);    
          
// c my correct end
//           printf("j %10d ub(j) %15.5e uf(j) %15.5e gb0 %15.5e vb0 %15.5e vf0 %15.5e pinv  %15.5e termx %15.5e\n",
// 		  j,     ub[j-1],     uf[j-1],     gb0,       vb0,       vf01,      pinv1,       termx);
	  
#ifdef DEBUG_INITIAL_PARTICLE_PRINTS
	       printf("electron %10d %25.15e %25.15e %25.15e %25.15e \n",2*j-2,yi[j - 1],uf[2*j-2],vf[2*j-2],wf[2*j-2]);
	       printf("electron %10d %25.15e %25.15e %25.15e %25.15e \n",2*j-1,yi[j - 1],uf[2*j-1],vf[2*j-1],wf[2*j-1]);
#endif
	  }

    free(ux);
    free(uy);
    free(uz);

#ifdef DEBUG_INITIAL_PARTICLE_PRINTS
//    exit(0);
#endif

    return 0;
} /* start_ */

int AddBeamParticles(int jmb,
				   double tex0,double tey0,double tez0,
				   double beam_lx, double beam_ly,int *jmb_real,
                   double lx,double ly,double lz, int meh,double Tb,double rimp,double rbd,
				   double *xb,double *yb, double *zb,double *ub,double *vb, double *wb
				  )
{
    double x,y,z,vb0,d__1,d__2,d__3,vy,vz,termx,gb0;
    double vf01,vf02,pinv1,pinv2,mfrq;
//     double *ux,*uy,*uz;
    double *ux,*uy,*uz;
    double beam_y_max,beam_y_min, beam_sh;

    beam_sh = (ly - beam_ly)/2;
    beam_y_max = ly - beam_sh;
    beam_y_min = beam_sh;

    for (int j = 1; j <= jmb; j++)
    {
    	z =   beam_ly * rnd_uniform(0) + beam_y_min;
	    y =   meh * ly + beam_ly * rnd_uniform(0) + beam_y_min;
	    x =   beam_lx * rnd_uniform(0);

	    xb[j - 1] = x;
	    yb[j - 1] = y;
	    zb[j - 1] = z;

	    double vb0,ux,uy,uz;
        vb0       = rnd_gaussian(0.0, Tb*rimp,0);
        ux        = vb0 + rimp;
        uy        = rnd_gaussian(0.0, Tb*rimp,0);
        uz        = rnd_gaussian(0.0, Tb*rimp,0);

        vb0 = sqrt(1.0 - ux*ux - uy * uy - uz * uz);

        ub[j - 1] = ux / vb0;
        vb[j - 1] = uy / vb0;
        wb[j - 1] = uz / vb0;


#ifdef ADD_BEAM_INITIAL_PARTICLE_PRINTS
	       printf("add beam %10d %25.15e  %25.15e    %25.15e %25.15e %25.15e   %25.15e \n",
	    		         j,  xb[j - 1],yb[j - 1],vb0,    ub[j-1],vb[j - 1],wb[j - 1]);
#endif
    }
	return 0;
}
/*     MAKING THE RANDOM GENERATOR WORK THE SAME */
//    parasetrandomelectrons_();
/*     END RANDOM GENERATOR */

//     vy = rnd_gaussian(0.0,tey0,1);
//     vz = rnd_gaussian(0.0,tez0,1);
//
//     termx = rnd_gaussian(0.0,tex0,1);
//    for (j = 1; j <= jmb;j++)
//    {
//           xf[2*j-1-1] = xi[j-1];
//           yf[2*j-1-1] = yi[j-1];
//           zf[2*j-1-1] = zi[j-1];
//
//           xf[2*j-1]   = xi[j-1];
//           yf[2*j-1]   = yi[j-1];
//           zf[2*j-1]   = zi[j-1];
//
////         FIRST SETTING TRANVERSE
//// razbros v skorostyax
//
//           vy=rnd_gaussian(0.0,tey0,0);
//           vz=rnd_gaussian(0.0,tez0,0);
//
////          INVERSE CURRENT
//
//           termx = rnd_gaussian(0.0,tex0,0);
//// 	   printf("termx %15.5e vx %15.5e vy %15.5e \n",termx,vy,vz);
//
//           gb0 = pow(1.0+pow(ub[j-1],2)+pow(vb[j-1],2)+pow(wb[j-1],2),-0.5);
//
//           vb0=ub[j-1]*gb0;
//      	   if ((beamf == 1) && ((xi[j - 1] < beam_lx) && (yi[j - 1] < beam_y_max) && (yi[j - 1] > beam_y_min)))
//     	   {
//        	    vf01=-rbd*vb0+termx;
//                vf02=-rbd*vb0-termx;
//	       }
//       	   else
//	       {
//      	        vf01=+termx;
//      	        vf02=-termx;
//	       }
//
//           pinv1= vf01*pow((1.0-pow(vf01,2)-vy*vy-vz*vz),-0.5);
//           pinv2= vf02*pow((1.0-pow(vf02,2)-vy*vy-vz*vz),-0.5);
//
//          vf[2*j-2] =  vy*pow((1.0-pow(vf01,2)-vy*vy-vz*vz),-0.5);
//// 	  write(37,127) 2*j-1,vy,vf01,vz,
////      +    vy/dsqrt((1d0-vf01**2-vy**2-vz**2))
////           printf("%10d vy %15.5e vf01 %15.5e vz %15.5e vf %15.5e \n",2*j-1,vy,vf01,vz,vy*pow((1.0-pow(vf02,2)-vy*vy-vz*vz),-0.5));
//          vf[2*j-1] = -vy*pow((1.0-pow(vf02,2)-vy*vy-vz*vz),-0.5);
//
//          wf[2*j-2] =  vz*pow((1.0-pow(vf01,2)-vy*vy-vz*vz),-0.5);
//          wf[2*j-1] = -vz*pow((1.0-pow(vf02,2)-vy*vy-vz*vz),-0.5);
//
//	      uf[2*j-2] = pinv1+0.01*sin(mfrq*2.0*M_PI*xf[2*j-2]/lx);
//          uf[2*j-1] = pinv2+0.01*sin(mfrq*2.0*M_PI*xf[2*j-1]/lx);

// c my correct end
//           printf("j %10d ub(j) %15.5e uf(j) %15.5e gb0 %15.5e vb0 %15.5e vf0 %15.5e pinv  %15.5e termx %15.5e\n",
// 		  j,     ub[j-1],     uf[j-1],     gb0,       vb0,       vf01,      pinv1,       termx);

//#ifdef DEBUG_INITIAL_PARTICLE_PRINTS
//	       printf("electron %10d %25.15e %25.15e %25.15e %25.15e \n",2*j-2,yi[j - 1],uf[2*j-2],vf[2*j-2],wf[2*j-2]);
//	       printf("electron %10d %25.15e %25.15e %25.15e %25.15e \n",2*j-1,yi[j - 1],uf[2*j-1],vf[2*j-1],wf[2*j-1]);
//#endif
//	  }

//    free(ux);
//    free(uy);
//    free(uz);
//
//#ifdef DEBUG_INITIAL_PARTICLE_PRINTS
////    exit(0);
//#endif
//
//    return 0;
//} /* start_ */




//double compare(double *a,double *b,int num,char *legend,double tol)
//{
//     double t = 0.0;
//     int    i;
//
//     for(i = 0; i < num ;i++)
//     {
//         if(fabs(a[i] - b[i]) < tol)
//         {
//            t += 1.0;
//#ifdef COMPARE_PRINTS
//            printf(" i %5d a %e b %e diff %e\n",i,a[i],b[i],fabs(a[i] - b[i]));
//#endif
//
//         }
//         else
//         {
//#ifdef COMPARE_PRINTS
//        	printf("WRONG i %5d a %e b %e diff %e\n",i,a[i],b[i],fabs(a[i] - b[i]));
//#endif
//         }
//     }
//
//     if(num > 0) t /= num;
//#ifdef COMPARE_PRINTS
//     printf("%30s %.5f\n",legend,t);
//#endif
//     return t;
//}

int getMassCharge(ParticleArrays *ions,ParticleArrays *electrons,ParticleArrays *beam_electrons,
		double ni,double rbd,int lp)
{
    //int lp = ((double)N)/(Nx*Ny*Nz);
	electrons->m[0]      = -ni/lp;                 //!!!!!!
	ions->m[0]           =  2.0*ni/lp*(1.0+rbd);
	beam_electrons->m[0] = electrons->m[0]*rbd*2.0;

	electrons->q_m        = -1.0;
	ions->q_m             =  1.0/1836.0;
	beam_electrons->q_m   = -1.0;

}


//int main(void)
//{
//  double tex0 = 1e-3,tey0 = 1e-3,tez0 = 1e-3;
//  double lx=1.1424,ly=0.05,lz=0.05,Tb = 0.14,rimp = 0.2,rbd = 2.0e-3;
//  int meh = 0;
//
//  InitUniformMaxwellianParticles(1,N,tex0,tey0,tez0,lx,ly,lz,meh,Tb,rimp,rbd,
//				   xi,yi, zi,ui,ui, wi,
//				   xb,yb, zb,ub,vb, wb,
//				   xf,yf, zf,uf,vf, wf
//				  );
//
//  char part_name[100];
//  int part_nt = 1;
//
//  sprintf(part_name,"mumu000%08d.dat",part_nt);
//
//
//  ParticleArrays ions,electrons,beam_electrons;
//
//  InitBinaryParticlesArrays(part_name,part_nt,&ions,&electrons,&beam_electrons,100,4,4);
//
//  double tol = 1e-15;
//  double t_ion_x = compare(xi,ions.dbg_x,N,"ionx",tol);
//  double t_ion_y = compare(yi,ions.dbg_y,N,"iony",tol);
//  double t_ion_z = compare(zi,ions.dbg_z,N,"ionz",tol);
//  double t_ion_px = compare(ui,ions.dbg_px,N,"ionpx",tol);
//  double t_ion_py = compare(vi,ions.dbg_py,N,"ionpy",tol);
//  double t_ion_pz = compare(wi,ions.dbg_pz,N,"ionpz",tol);
//  printf("IONS x %15.6e y %15.6e z %15.6e px %15.6e py %15.6e pz %15.6e \n",t_ion_x,t_ion_y,t_ion_z,t_ion_px,t_ion_py,t_ion_pz);
//
//  double t_electron_x = compare(xf,electrons.dbg_x,2*N,"electronx",tol);
//  double t_electron_y = compare(yf,electrons.dbg_y,2*N,"electrony",tol);
//  double t_electron_z = compare(zf,electrons.dbg_z,2*N,"electronz",tol);
//  double t_electron_px = compare(uf,electrons.dbg_px,2*N,"electronpx",tol);
//  double t_electron_py = compare(vf,electrons.dbg_py,2*N,"electronpy",tol);
//  double t_electron_pz = compare(wf,electrons.dbg_pz,2*N,"electronpz",tol);
//  printf("ELECTRONS x %15.6e y %15.6e z %15.6e px %15.6e py %15.6e pz %15.6e \n",t_electron_x,t_electron_y,t_electron_z,t_electron_px,t_electron_py,t_electron_pz);
//
//  double t_beam_electron_x = compare(xb,beam_electrons.dbg_x,N,"beam_electronx",tol);
//  double t_beam_electron_y = compare(yb,beam_electrons.dbg_y,N,"beam_electrony",tol);
//  double t_beam_electron_z = compare(zb,beam_electrons.dbg_z,N,"beam_electronz",tol);
//  double t_beam_electron_px = compare(ub,beam_electrons.dbg_px,N,"beam_electronpx",tol);
//  double t_beam_electron_py = compare(vb,beam_electrons.dbg_py,N,"beam_electronpy",tol);
//  double t_beam_electron_pz = compare(wb,beam_electrons.dbg_pz,N,"beam_electronpz",tol);
//  printf("BEAM x %15.6e y %15.6e z %15.6e px %15.6e py %15.6e pz %15.6e \n",t_beam_electron_x,t_beam_electron_y,t_beam_electron_z,t_beam_electron_px,
//	 t_beam_electron_py,t_beam_electron_pz);
//
//
//  return 0;
//}
