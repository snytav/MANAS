/*
 * xeon_plasma.h
 *
 *  Created on: Oct 30, 2015
 *      Author: snytav
 */

#ifndef XEON_PLASMA_H_
#define XEON_PLASMA_H_

#include "gpu_plasma.h"


template <template <class Particle,int dims> class Cell,int dims >
class XeonPlasma: public GPUPlasma<Cell,dims>
{
public:
	 XeonPlasma(int nx,int ny,int nz,double lx,double ly,double lz,
			 int n_per_cell1,double q_m,double TAU,int f3D,double Bx0,int N0,
			 double tex0,double tey0,double tez0,double Tb,double rimp,double rbd,double ni,
			 double bx,double by,double bz,
			 int beam_plasma,int start_from_file):
		     GPUPlasma<Cell,dims>(nx,ny,nz,lx,ly,lz,n_per_cell1,q_m,TAU,f3D,Bx0,N0,
				 tex0,tey0,tez0,Tb,rimp,rbd,ni,
				 bx,by,bz,beam_plasma,start_from_file){}

	 ~XeonPlasma(){}
	 XeonPlasma(){}

};

#endif /* XEON_PLASMA_H_ */
