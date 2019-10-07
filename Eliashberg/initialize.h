#ifndef _initialize_H
#define _initialize_H

#include "some_function.h"

complex_3D initialize_complex()
{
	int ix,iy,ie;
	complex_3D C ;
	C.s1=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.x1=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.p1=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.s2=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.x2=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.p2=d3tensor2(0,Nener,0,nmom,0,nmom);
	for(ie=0;ie<=Nener;ie++)
	{
		for(ix=0;ix<=nmom;ix++)
		{
			for(iy=0;iy<=ix;iy++)
			{
				C.s1[ie][ix][iy] = 0;
                C.x1[ie][ix][iy] = 0;
                C.p1[ie][ix][iy] = 0;
                C.s2[ie][ix][iy] = 0;
                C.x2[ie][ix][iy] = 0;
                C.p2[ie][ix][iy] = 0;
			}
		}
	}
	return C;
}
void free_complex_3D(complex_3D C)
{
	free_d3tensor(C.s1,0,Nener,0,nmom,0,nmom);
	free_d3tensor(C.s2,0,Nener,0,nmom,0,nmom);
	free_d3tensor(C.x1,0,Nener,0,nmom,0,nmom);
	free_d3tensor(C.x2,0,Nener,0,nmom,0,nmom);
	free_d3tensor(C.p1,0,Nener,0,nmom,0,nmom);
	free_d3tensor(C.p2,0,Nener,0,nmom,0,nmom);

}
void free_double_3D(double_3D A)
{
	free_d3tensor(A.s,0,Nener,0,nmom,0,nmom);
	free_d3tensor(A.x,0,Nener,0,nmom,0,nmom);
	free_d3tensor(A.p,0,Nener,0,nmom,0,nmom);

}

complex_3D initialize_eff()
{
	int ix,iy,ie;
	double kx, ky;
	complex_3D C;
	C.s1=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.x1=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.p1=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.s2=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.x2=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.p2=d3tensor2(0,Nener,0,nmom,0,nmom);
	for(ie=0;ie<=Nener;ie++)
	{
		for(ix=0;ix<=nmom;ix++)
		{
			kx = momentum(ix);
			for(iy=0;iy<=ix;iy++)
			{
				ky = momentum(iy);
				C.s1[ie][ix][iy] = 0;
                C.x1[ie][ix][iy] = 0;
                C.p1[ie][ix][iy] = 0.015*d_wave(kx,ky);
                C.s2[ie][ix][iy] = 0;
                C.x2[ie][ix][iy] = 0;
                C.p2[ie][ix][iy] = 0;
			}
		}
	}
	return C;
}

complex_3D initialize_imp()
{
	int ix,iy,ie;
	complex_3D C;
	C.s1=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.x1=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.p1=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.s2=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.x2=d3tensor2(0,Nener,0,nmom,0,nmom);
	C.p2=d3tensor2(0,Nener,0,nmom,0,nmom);
	for(ie=0;ie<=Nener;ie++)
	{
		for(ix=0;ix<=nmom;ix++)
		{
			for(iy=0;iy<=ix;iy++)
			{
				C.s1[ie][ix][iy] = 0;
                C.x1[ie][ix][iy] = 0;
                C.p1[ie][ix][iy] = 0;
                C.s2[ie][ix][iy] = -small*t1;
                C.x2[ie][ix][iy] = 0;
                C.p2[ie][ix][iy] = 0;
			}
		}
	}
	return C;
}



double_3D initialize_real()
{
	int ix,iy,ie;
	double_3D A;
	A.s=d3tensor2(0,Nener,0,nmom,0,nmom);
	A.x=d3tensor2(0,Nener,0,nmom,0,nmom);
	A.p=d3tensor2(0,Nener,0,nmom,0,nmom);
	for(ie=0;ie<=Nener;ie++)
	{
		for(ix=0;ix<=nmom;ix++)
		{
			for(iy=0;iy<=ix;iy++)
			{
				A.s[ie][ix][iy] = 0;
                A.x[ie][ix][iy] = 0;
                A.p[ie][ix][iy] = 0;
			}
		}
	}
	return A;
}

double **initialize_Vq_imp()
{

	int ix, iy, ix2,iy2;
	double kx, ky;
	double **Vq_imp;

	Vq_imp = dmatrix(0,nmom,0,nmom);
	for(ix=0;ix<=nmom;ix++)
	{
		kx = momentum(ix)  ;
		for(iy=0;iy<=ix;iy++)	
		{
			ky = momentum(iy) ;
			Vq_imp[ix][iy] = niv0*pow(t1,2)*pow(Vimp(kx,ky),2)*pow(kapa_imp,3) ;
		}
	}

	fourn_m_to_r(Vq_imp,1);

	return Vq_imp;
}



double ***initialize_Vq_spin()
{

	FILE *f_delta, *f_kapa, *f_xamp;
	int ix, iy, ix2,iy2,iwb;
	double kx, ky, wboson ;
	double *delta, *kapa, *xamp;
	double lambda, Rq ;
	double ***Vq_eff;


	delta = dvector(0,Nener);
	kapa = dvector(0,Nener);
	xamp = dvector(0,Nener);

	for(iwb=0;iwb<=Nener;iwb++)
	{
		delta[iwb] = 0.0 ;
		kapa[iwb]= 0.0 ;
		xamp[iwb]= 0.0 ;
	}
	
	f_delta = fopen("delta.d","r");
	f_kapa = fopen("kapa.d","r");
	f_xamp = fopen("xamp.d","r");

	for(iwb=1;iwb<=iwbc;iwb++)
	{
		fscanf(f_delta,"%lf\t%lf\n",&wboson,&delta[iwb]);
		fscanf(f_kapa,"%lf\t%lf\n",&wboson,&kapa[iwb]);
		fscanf(f_xamp,"%lf\t%lf\n",&wboson,&xamp[iwb]);
	}
	fclose(f_delta) ;
	fclose(f_kapa) ;
	fclose(f_xamp) ;

	Vq_eff = d3tensor2(0,iwbc,0,nmom,0,nmom);
	for(iwb=1;iwb<=iwbc;iwb++)
	{
		wboson = energy_boson(iwb) ;
		if(wboson>0.040)
		{
			lambda = 0 ;
		}
		else
		{
			lambda = 4 ;
		}
		for(ix=0;ix<=nmom;ix++)
		{
			kx = momentum(ix)  ;
			for(iy=0;iy<=ix;iy++)	
			{
				ky = momentum(iy)   ;
				
				Rq = (pow(pow(kx-1.0,2)+pow(ky-1.0,2)-pow(2*delta[iwb],2),2)+ lambda*pow(kx-1.0,2)*pow(ky-1.0,2))/(4*pow(4*delta[iwb],2)) ;
				
				Vq_eff[iwb][ix][iy] = intensity*0.005*xamp[iwb]* pow(kapa[iwb],4)/pow(pow(kapa[iwb],2)+Rq,2);
			}
		}
	}

	free_dvector(delta,0,Nener);
	free_dvector(kapa,0,Nener);
	free_dvector(xamp,0,Nener);


	for(iwb=1;iwb<=iwbc;iwb++)
	{
		fourn_m_to_r(Vq_eff[iwb],1);
	}

	return Vq_eff;
}

#endif

