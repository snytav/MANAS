/*
 * control.h
 *
 *  Created on: Apr 26, 2014
 *      Author: snytav
 */

#ifndef CONTROL_H_
#define CONTROL_H_

#include <stdio.h>



#define FLY_PRINTS

//typedef enum {ION, PLASMA_ELECTRON, BEAM_ELECTRON} particle_sorts;
#include "particle_types.h"

#ifdef __CUDACC__
__host__ __device__
#endif
void writeParticleAttribute(int j,double ami,int num,double t);

#ifdef __CUDACC__
__host__ __device__
#endif
int ParticleAttributePosition(int jmp,int j,int sort,int num)
{
    int pos;//,pos1,pos2,pos0;

    if(sort == ION)
    {
       pos = num+(j-1)*PARTICLE_ATTRIBUTES ;
    }

    if(sort == PLASMA_ELECTRON)
    {
       pos = num+(j-1)*PARTICLE_ATTRIBUTES + jmp*PARTICLE_ATTRIBUTES ;
//       pos1 = pos;
    }

    if(sort == BEAM_ELECTRON)
    {
       pos = num+(j-1)*PARTICLE_ATTRIBUTES + jmp*2*PARTICLE_ATTRIBUTES ;
//       pos2 = pos;
    }

//    pos0 = num+(j-1)*PARTICLE_ATTRIBUTES ;
//    pos1 = num+(j-1)*PARTICLE_ATTRIBUTES + jmp*PARTICLE_ATTRIBUTES ;
//    pos2 = num+(j-1)*PARTICLE_ATTRIBUTES + jmp*2*PARTICLE_ATTRIBUTES ;
//
//    //printf("PAP %10d %10d %10d %10d %10d %10d \n",jmp,j,num,pos0,pos1,pos2);

    //printf("sort %d pos %d \n",(int)sort,pos);
    return pos;
}

#ifdef __CUDACC__
__host__ __device__
#endif
int ParticleAttributePositionFortran(int jmp,int j,int sort,int num)
{

	return ParticleAttributePosition(jmp,j,sort,num) - 1;
}


#endif /* CONTROL_H_ */
