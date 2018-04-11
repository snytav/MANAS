/*
 * diagnose.c
 *
 *  Created on: Jun 5, 2016
 *      Author: snytav
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "plasma/physics.h"
#include "diagnose.h"
#include "particle_types.h"

int flag_output_allparticles = 1;


int output_kernel(float *s2,float *s4,int *i,int *l,float x,float y,float hx,float hy)
{
	*s2 = x/hx;
	*i  = (int)(*s2+1.50);
	*s2 = *s2+1.50-*i;
	*s4 = y/hy;
	*l  = (int)(*s4+1.50);
	*s4 = *s4+1.50-*l;
}

float particle_physical_characteristic(
		float m,
		float q_m,
		float pu,
		float pv,
		float pw,
		particle_diagnostic_operation code
		)
{
    double ps=1.0/sqrt(1.0+(pu*pu+pv*pv+pw*pw));
    double ms = m,a;

    switch(code)
    {
    case CHARGE: a = ms;
                 break;

    case VELX:   a = ps*pu;
                 break;

    case VELY:   a = ps*pv;
                 break;

    case VELZ:   a = ps*pw;
                 break;
    }

    return a;
}

int write_to_matrix(diagnose_matrix m,float a,int i,int l,float s2,float s4)
{
	float s,s21,s41;

	i--;
	l--;

	s21 = 1.0-s2;
	s41 = 1.0-s4;

    s = s21*s41*a;
    m[i][l] += s;

    s=s21*s4*a;
    m[i][l+1] += s;

    s=s2*s41*a;
    m[i+1][l] += s;

    s=s2*s4*a;
    m[i+1][l+1] += s;

    return 0;
}

int write_matrix_to_file(char *shot_name,int nt,char *name,
		                 diagnose_matrix m,float hx,float hy)
{
	char fname[200];
	FILE *f;
	int i,l;


	sprintf(fname,"%s_nt%08d_%s.dat",name,nt,shot_name);
	if((f = fopen(fname,"wt")) == NULL) return 1;

	for(i = 0;i <DIAGNOSE_NX;i++)
	{
		for(l = 0;l < DIAGNOSE_NY;l++)
		{
			fprintf(f,"%15.5f %15.5f %15.5e \n",i*hx,l*hy,m[i][l]);
		}
		fprintf(f,"\n");
	}
	fclose(f);

	return 0;
}

int TransformParticlesToLists(int sorts,
        ParticleArraysGroup pag,
        float *particles,
        ParticleFloatArraysGroup pfag
		)
{
	int is,fn,total;

	for(is = 0;is < sorts;is++)
	{
		total = (pag[is]).total;
		for(fn = 0;fn < pag[is].total;fn++)
		{
            pfag[is].dbg_x[fn]   = particles[8*total*is+fn*8    ];
            pfag[is].dbg_y[fn]   = particles[8*total*is+fn*8 + 1];
            pfag[is].dbg_z[fn]   = particles[8*total*is+fn*8 + 2];
            pfag[is].dbg_px[fn]  = particles[8*total*is+fn*8 + 3];
            pfag[is].dbg_py[fn]  = particles[8*total*is+fn*8 + 4];
            pfag[is].dbg_pz[fn]  = particles[8*total*is+fn*8 + 5];
		}
		pfag[is].m     = pag[is].m[0];
		pfag[is].q_m   = pag[is].q_m;
		pfag[is].total = pag[is].total;
	}

	return 0;
}


int Particle2D_Distributions(int sorts,char *shot_name,int nt,
		 ParticleFloatArraysGroup pfag,
		 double Lx,double Ly
		 )
{
	float hx = Lx/DIAGNOSE_NX,hy = Ly/DIAGNOSE_NY;
	float x,y,pu,pv,pw,m,q_m,z;
	float n[DIAGNOSE_NX][DIAGNOSE_NY];
	int i,k,is,fn;

	memset(n,0,sizeof(float)*DIAGNOSE_NX*DIAGNOSE_NY);

	for(is = 0;is < sorts;is++)
	{
		m   = pfag[is].m;
		q_m = pfag[is].q_m;

		for(fn = 0;fn < pfag[is].total;fn++)
		{
		    float s2,s4,a;
            int i,l;

            x  = pfag[is].dbg_x[fn];
            y  = pfag[is].dbg_y[fn];
            pu = pfag[is].dbg_px[fn];
            pv = pfag[is].dbg_py[fn];
            pw = pfag[is].dbg_pz[fn];

            output_kernel(&s2,&s4,&i,&l,x,y,hx,hy);

            a =  particle_physical_characteristic(m,q_m,pu,pv,pw,CHARGE);

            write_to_matrix(n,a,i,l,s2,s4);
		}
	}
//	write_matrix_to_file(shot_name,nt,"n",n,hx,hy);



	return 0;
}

int write_particles_component(
		float *x,float *y,float *z,float *pu,float *pv, float *pw,
		int num, char *comp_name, char *shot_name, int nt
		)
{
	int j;
	char fname[1000];
	FILE *f;

	sprintf(fname,"%s_%s_%08d.dat",comp_name,shot_name,nt);

	if((f = fopen(fname,"wt")) == NULL) return 1;

	for(j = 0;j < num;j++)
	{
		fprintf(f,"%10d %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e\n",j,
				x[j],
				y[j],
				z[j],
				pu[j],
				pv[j],
				pw[j]
				);
	}
	fclose(f);

}

int write_particles(ParticleFloatArraysGroup *pfag,ParticleArraysGroup *init,char *variant_name,int nt)
{
	int i;
	if(flag_output_allparticles != 1) return 0;

	for(i = 0;i < SORTS;i++)
	{
		char sort_name[100];

		particle_sort_name(sort_name,(particle_sorts)i);

	    write_particles_component(
			 		(*pfag)[i].dbg_x,(*pfag)[i].dbg_y,(*pfag)[i].dbg_z,
			 		(*pfag)[i].dbg_px,(*pfag)[i].dbg_py, (*pfag)[i].dbg_pz,
			 		(*pfag)[i].total,sort_name, variant_name,nt);
	}

    return 0;
}

int particle_sort_name(char *name,particle_sorts sort)
{
	switch(sort)
	{
	case ION:             strcpy(name,"ion");
 	                      break;

	case PLASMA_ELECTRON: strcpy(name,"electron");
 	                      break;

	case BEAM_ELECTRON:   strcpy(name,"beam");
 	                      break;
	}
}



