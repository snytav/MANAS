/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdio.h>
#include <stdlib.h>



int readParameterFile(int *beam_plasma,int *start_from_file,
		              double *tex0,double *tey0,double *tez0,
		              double *Tb,double *rimp,
		              double *rbd,double *ni,
			          double *lx,double *ly,double *lz,
			          int *lp,int *nx,int *ny,int *nz,
			          double *tau,double *B0,
			          int *np)
{
	FILE *f;
	char str[1000];
	int n;

	if((f = fopen("0_plasma.dat","rt")) == NULL)
	{
		return 1;
	}
	fgets(str,1000,f);
	*beam_plasma = atoi(str);

	fclose(f);


	if(beam_plasma == 0)
	{
		chdir("./beam");
	}
	else
	{
		chdir("./plasma");
	}

	if((f = fopen("000_params.dat","rt")) == NULL)
	{
		return 1;
	}



	return 0;
}
