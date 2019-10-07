#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "nrutil.c"
#include "fourn.c"

/* optimal
#define t1 0.420
#define t2 0.073
#define t3 0.036
#define chem -0.3200
*/
//over
#define t1 0.250
#define t2 0.0375
#define t3 0.01875
#define chem -0.190
#define tol 0.00000001
#define tol2 0.00000001
#define tol3 0.00000001
#define coupling 0.2
#define gamma 0.15
#define gamma2 0.025 

double niv0 = 0.03 ;
double kapa_imp = 0.3 ;
double V0 = 1.0 ;

//double intensity = 0.560346*pow(2*M_PI,2) ;
double intensity = 22.1216 ;

double Temp = 0.00008617*000.0 ;
//double Temp = 0.000000001 ;
double wb = 0.050 ;
int nmom = 128, nmom5 = 500, Nener = 3000 ;
double enstep, small, step, step2 ;
double wc = 3 ;
int iwb, iwbs = 1, iwbc = 400  ;

double xi_band(double kx, double ky) ;
double fermi_d(double w) ;
double bose_d(double wboson) ;
double d_wave(double kx, double ky) ;
//double asq_F(double wboson,double kx, double ky) ;
double Vimp(double qx, double qy);


int main(void)
{
	FILE *fp,*fp0,*fp00, *fp1, *fp111, *fp2, *fp200, *fp22, *fp202,*fp23,*fp24,*fp25,*fp26,*fp27,*fp28, *fp3, *fp4, *fp999, *fp6, *fp7, *fp11,*fp12,*fp30, *fp31, *fp35, *fp50,*fp_dre,*fp_dre1,*fp_dre2,*fp_dim, *fp_p, *fp8;
	FILE *f_delta, *f_kapa, *f_xamp, *fp61, *fp62, *fp63, *fp64, *fpp, *fpp0, *fpp5, *fpp10, *fpp15, *fpp20, *fpp25, *fpp30, *fpp35, *fpp40,*fp_z,*fp_d,*fp_df,*fpp00,*fp_edc,*fp_angle,*fp_angle2, *fp_angle3 ; 
	FILE *fpmi, *fp_sep, *fp_sep_s, *fp_sep_x, *fp_sep_p ;
	char fid[20],fid2[20],fid3[20],fid4[20],fid5[20],fid000[20],fid_mdc[20],fid_edc[20],fid_angle[20],fid_angle2[20],fid_angle3[20] ;	
	char fi_se[40],fi_se2[40],fi_se3[40],fi_se4[40],fi_se5[40],fi_se000[40],fi_mdc[40],fi_edc[40],fi_angle[40],fi_angle2[40],fi_angle3[40] ;

	unsigned long *nmom2 ;
	int i, j, k, l , ie, ie0, iee, je, ien, ied, ix, iy, ix1,iy1, ix2, iy2, ix3,iy3, ix0, iy0, iter, iter2, iter22, iter3, N_iter, N_iter2, N_iter22, N_iter3, re_temp , swi, iek ;
	int ishift , itheta, theta ;
	double x, y, kx, ky, kx1,ky1, kx2,ky2, kx3, ky3, kx0,ky0, w, w2, error, error_s, error_x, error_p , xi_kxky, ener, wboson, fermi, bose, min_band , ek , e_temp, doping,dop_temp, dek, det, aa ;
	double numer1, numer2, deno1, deno2, err2, error2, error3, err3,err_s, err_x, err_p , delt, dei, enk,den_1,den_2 ,den0, R, R2 , cpi , peak, comp ;
	double  pair_re, pair_im, imp_s1, imp_s2, imp_x1, imp_x2, imp_p1, imp_p2 ;
	float *ws1, *xx1, *nu_r, *nu_i, *de_r, *de_i ;
	double deltkk, Rq, lamda, Amax ;
	double temp1, temp2 ;
	double d_re_in, d_re_co;
	double count_x,count_y, temp, mixing, mixing3;
	double A11,A13,A31,A33, s1kxky, s2kxky, x1kxky, x2kxky, p1kxky, p2kxky;
	int Ami ;

	float ***self_imp_re, ***self_imp_im, ***self_imp_re_new, ***self_imp_im_new;	
	float ***renor_imp_re, ***renor_imp_im, ***renor_imp_re_new, ***renor_imp_im_new ;
	float ***pairing_imp_re, ***pairing_imp_im, ***pairing_imp_re_new, ***pairing_imp_im_new;


	int *band_0, *band_delta ;
	float ***den,***den_new , *dos, *dos_im, *dos_im2;
	double *delta, *kapa, *xamp , **z_re_0, **d_re_0, *z_re_e, *z_im_e, *d_re_e, **Akxky ;
	float ***self_re, ***self_im, ***self_re_new, ***self_im_new, *self_temp ;	
	float ***renor_re, ***renor_im, ***renor_re_new, ***renor_im_new, *renor_temp ;
	float ***pairing_re, ***pairing_im, ***pairing_re_new, ***pairing_im_new, *pairing_temp ;
	float *s_new_re, *s_new_im, *x_new_re, *x_new_im, *p_new_re, *p_new_im , *pair_re_e, *pair_im_e ;
	float *asqf, *Akxky0 ;
	float ***pairing_re_in, ***pairing_im_in, ***pairing_re_co, ***pairing_im_co ;

	float ***Gs_re, ***Gx_re, ***Gp_re, ***Gs_im, ***Gx_im, ***Gp_im ;
	float ***Vq, *Vq_temp;
	float *As_temp, ***As, *Ax_temp, ***Ax, *Ap_temp, ***Ap;
	float *Ns_re, *Ns_im, *Nx_re, *Nx_im ; 

	float **Vq_imp, *Vq_imp_temp;
	float *Gs_temp, *Gx_temp, *Gp_temp ;
	
	float ***sep_s_re,***sep_s_im,***sep_x_re,***sep_x_im,***sep_p_re,***sep_p_im ;
	float *s_sep_re,*s_sep_im, *x_sep_re, *x_sep_im, *p_sep_re, *p_sep_im ;
	float *sep_s_temp, *sep_x_temp, *sep_p_temp ;

	time_t time_0, time_now, running_time ;



	fp = fopen("mat.d","w") ;
// time
	time(&time_0);
	printf(ctime(&time_0));
	fprintf(fp,ctime(&time_0));

	N_iter = 15 ;
	N_iter2 = 10 ;
	N_iter22 = 0 ;
	N_iter3 = 10 ;


	enstep = (1.0*wc)/Nener ;
	small = 3.0*enstep/1.0 ;
	step = 1.0/nmom ;
	step2 = 1.0/nmom5 ;
	delt = -3.0*enstep/1.0 ;
	deltkk = enstep/10000.0 ;
//	wb = iwb * enstep ;
	printf("%f\t%f\n",enstep , small) ;

	nmom2 = lvector(1,2) ;
	nmom2[1] = 2*nmom ;
	nmom2[2] = 2*nmom ;

	self_re = f3tensor(0,Nener,0,nmom,0,nmom) ;
	self_im = f3tensor(0,Nener,0,nmom,0,nmom) ;
	self_re_new = f3tensor(0,Nener,0,nmom,0,nmom) ;
	self_im_new = f3tensor(0,Nener,0,nmom,0,nmom) ;
	renor_re = f3tensor(0,Nener,0,nmom,0,nmom) ;
	renor_im = f3tensor(0,Nener,0,nmom,0,nmom) ;
	renor_re_new = f3tensor(0,Nener,0,nmom,0,nmom) ;
	renor_im_new = f3tensor(0,Nener,0,nmom,0,nmom) ;
	pairing_re = f3tensor(0,Nener,0,nmom,0,nmom) ;
	pairing_im = f3tensor(0,Nener,0,nmom,0,nmom) ;
	pairing_re_new = f3tensor(0,Nener,0,nmom,0,nmom) ;
	pairing_im_new = f3tensor(0,Nener,0,nmom,0,nmom) ;
	pairing_re_in = f3tensor(0,Nener,0,nmom,0,nmom) ;
	pairing_im_in = f3tensor(0,Nener,0,nmom,0,nmom) ;
	pairing_re_co = f3tensor(0,Nener,0,nmom,0,nmom) ;
	pairing_im_co = f3tensor(0,Nener,0,nmom,0,nmom) ;

	self_temp = vector(1,2*4*nmom*nmom) ;
	renor_temp = vector(1,2*4*nmom*nmom) ;
	pairing_temp = vector(1,2*4*nmom*nmom) ;
	
//	den = f3tensor(0,Nener,1,nmom,1,nmom) ;
//	den_new = f3tensor(0,Nener,1,nmom,1,nmom) ;
	As = f3tensor(0,Nener,0,nmom,0,nmom) ;
	Ax = f3tensor(0,Nener,0,nmom,0,nmom) ;
	Ap = f3tensor(0,Nener,0,nmom,0,nmom) ;

	As_temp = vector(1,2*4*nmom*nmom) ;
	Ax_temp = vector(1,2*4*nmom*nmom) ;
	Ap_temp = vector(1,2*4*nmom*nmom) ;

	Vq = f3tensor(0,Nener,0,nmom,0,nmom) ;
	Vq_temp = vector(1,2*4*nmom*nmom) ; 
	Vq_imp = matrix(0,nmom,0,nmom) ;
	Vq_imp_temp = vector(1,2*4*nmom*nmom) ; 

	z_re_0 = dmatrix(0,nmom,0,nmom) ;
	d_re_0 = dmatrix(0,nmom,0,nmom) ;
	
	Akxky = dmatrix(1,6001,0,nmom5) ;
	Akxky0 = vector(0,nmom5) ;
	
	z_re_e = dvector(0,Nener) ;
	z_im_e = dvector(0,Nener) ;
	d_re_e = dvector(0,Nener) ;

	s_new_re = vector(0,Nener) ;
	s_new_im = vector(0,Nener) ;
	x_new_re = vector(0,Nener) ;
	x_new_im = vector(0,Nener) ;
	p_new_re = vector(0,Nener) ;
	p_new_im = vector(0,Nener) ;

	Ns_re = vector(0,Nener) ;
	Ns_im = vector(0,Nener) ;
	Nx_re = vector(0,Nener) ;
	Nx_im = vector(0,Nener) ;
		
	pair_re_e = vector(0,Nener) ;
	pair_im_e = vector(0,Nener) ;
	
	band_0 = ivector(0,nmom) ;
	band_delta = ivector(0,nmom) ;

	dos = vector(0,Nener) ;
	dos_im = vector(0,Nener) ;
	dos_im2 = vector(0,Nener) ;

	delta = dvector(0,Nener) ;
	kapa = dvector(0,Nener) ;
	xamp = dvector(0,Nener) ;

	self_imp_re = f3tensor(0,Nener,0,nmom,0,nmom) ;
	self_imp_im = f3tensor(0,Nener,0,nmom,0,nmom) ;
	self_imp_re_new = f3tensor(0,Nener,0,nmom,0,nmom) ;
	self_imp_im_new = f3tensor(0,Nener,0,nmom,0,nmom) ;
	renor_imp_re = f3tensor(0,Nener,0,nmom,0,nmom) ;
	renor_imp_im = f3tensor(0,Nener,0,nmom,0,nmom) ;
	renor_imp_re_new = f3tensor(0,Nener,0,nmom,0,nmom) ;
	renor_imp_im_new = f3tensor(0,Nener,0,nmom,0,nmom) ;
	pairing_imp_re = f3tensor(0,Nener,0,nmom,0,nmom) ;
	pairing_imp_im = f3tensor(0,Nener,0,nmom,0,nmom) ;
	pairing_imp_re_new = f3tensor(0,Nener,0,nmom,0,nmom) ;
	pairing_imp_im_new = f3tensor(0,Nener,0,nmom,0,nmom) ;

	Gs_re = f3tensor(0,Nener,0,nmom,0,nmom) ;
	Gx_re = f3tensor(0,Nener,0,nmom,0,nmom) ;
	Gp_re = f3tensor(0,Nener,0,nmom,0,nmom) ;
	Gs_im = f3tensor(0,Nener,0,nmom,0,nmom) ;
	Gx_im = f3tensor(0,Nener,0,nmom,0,nmom) ;
	Gp_im = f3tensor(0,Nener,0,nmom,0,nmom) ;

	Gs_temp = vector(1,2*4*nmom*nmom) ;
	Gx_temp = vector(1,2*4*nmom*nmom) ;
	Gp_temp = vector(1,2*4*nmom*nmom) ;

	ws1 = vector(0,Nener) ;
	xx1 = vector(0,Nener) ;
//	nu_r = vector(0,Nener) ;
//	nu_i = vector(0,Nener) ;
	de_r = vector(0,Nener) ;
	de_i = vector(0,Nener) ;

	sep_s_re = f3tensor(0,Nener,0,nmom,0,nmom) ;
	sep_s_im = f3tensor(0,Nener,0,nmom,0,nmom) ;
	sep_x_re = f3tensor(0,Nener,0,nmom,0,nmom) ;
	sep_x_im = f3tensor(0,Nener,0,nmom,0,nmom) ;
	sep_p_re = f3tensor(0,Nener,0,nmom,0,nmom) ;
	sep_p_im = f3tensor(0,Nener,0,nmom,0,nmom) ;

	s_sep_re = vector(0,Nener) ;
	s_sep_im = vector(0,Nener) ;
	x_sep_re = vector(0,Nener) ;
	x_sep_im = vector(0,Nener) ;
	p_sep_re = vector(0,Nener) ;
	p_sep_im = vector(0,Nener) ;

	sep_s_temp = vector(1,2*4*nmom*nmom) ;
	sep_x_temp = vector(1,2*4*nmom*nmom) ;
	sep_p_temp = vector(1,2*4*nmom*nmom) ;






	#pragma omp parallel for private(je)
	for(je=0;je<=Nener;je++) 
	{
		for(ix=0;ix<=nmom;ix++)
		{
			for(iy=0;iy<=nmom;iy++)	
			{
				self_imp_re[je][ix][iy] = 0 ;
				self_imp_im[je][ix][iy] = -small*t1 ;
				renor_imp_re[je][ix][iy] = 0 ;
				renor_imp_im[je][ix][iy] = 0 ;
				pairing_imp_re[je][ix][iy] = 0 ;
				pairing_imp_im[je][ix][iy] = 0 ;
				self_imp_re_new[je][ix][iy] = 0 ;
				self_imp_im_new[je][ix][iy] = 0 ;
				renor_imp_re_new[je][ix][iy] = 0 ;
				renor_imp_im_new[je][ix][iy] = 0 ;
				pairing_imp_re_new[je][ix][iy] = 0 ;
				pairing_imp_im_new[je][ix][iy] = 0 ;
			}
		}
	}

	for(ix=0;ix<=nmom;ix++)
	{
		kx = step*(ix) ;
		for(iy=0;iy<=nmom;iy++)	
		{
			ky = step*(iy) ;

			for(ie=0;ie<=Nener;ie++) 
			{

				self_re[ie][ix][iy] = 0.00 ;
				self_im[ie][ix][iy] = -0.00 ;
				renor_re[ie][ix][iy] = 0.00 ;
				renor_im[ie][ix][iy] = 0.00 ;
//				pairing_re[ie][ix][iy] = 0.00 ;
				pairing_re[ie][ix][iy] = 0.015*d_wave(kx,ky) ;
				pairing_im[ie][ix][iy] = 0.00 ;
				
			}
		}
	}


	for(ix=0;ix<=nmom;ix++)
	{
		kx = step*(ix)  ;
		for(iy=0;iy<=nmom;iy++)	
		{
			ky = step*(iy) ;

			Vq_imp[ix][iy] = niv0*pow(t1,2)*pow(Vimp(kx,ky),2)*pow(kapa_imp,3) ;

		}
	}

// fourier Vq_imp start
	for(ix=1;ix<=2*nmom;ix++)
	{
		for(iy=1;iy<=2*nmom;iy++)	
		{
			if(ix<=nmom)
			{
				ix2 = ix - 1 ;
			}
			else
			{
				ix2 = 2*nmom - ix + 1 ;
			}
			if(iy<=nmom)
			{
				iy2 = iy - 1 ;
			}
			else
			{
				iy2 = 2*nmom - iy + 1 ;
			}
			Vq_imp_temp[(2*ix-1)+(iy-1)*2*2*nmom]= Vq_imp[ix2][iy2] ; 
			Vq_imp_temp[2*ix+(iy-1)*2*2*nmom] = 0 ;
		}
	}
	fourn(Vq_imp_temp,nmom2,2,-1) ;
	for(ix=0;ix<=nmom;ix++)
	{
		for(iy=0;iy<=nmom;iy++)	
		{
			ix2 = ix + 1 ;
			iy2 = iy + 1 ;
			Vq_imp[ix][iy] = Vq_imp_temp[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
		}
	}



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

	for(iwb=1;iwb<=iwbc;iwb++)
	{
		wboson = iwb * enstep ;
		if(wboson>0.040)
		{
			lamda = 0 ;
		 }
		else
		{
			lamda = 4 ;
		}
		for(ix=0;ix<=nmom;ix++)
		{
			kx = step*(ix)  ;
			for(iy=0;iy<=nmom;iy++)	
			{
				ky = step*(iy)   ;
				
				Rq = (pow(pow(kx-1.0,2)+pow(ky-1.0,2)-pow(2*delta[iwb],2),2)+ lamda*pow(kx-1.0,2)*pow(ky-1.0,2))/(4*pow(4*delta[iwb],2)) ;
				
				Vq[iwb][ix][iy] = intensity*0.005*xamp[iwb]* pow(kapa[iwb],4)/pow(pow(kapa[iwb],2)+Rq,2);
			}
		}
	}



// fourier vq start
	for(iwb=1;iwb<=Nener;iwb++)
	{
		for(ix=1;ix<=2*nmom;ix++)
		{
			for(iy=1;iy<=2*nmom;iy++)	
			{
				if(ix<=nmom)
				{
					ix2 = ix - 1 ;
				}
				else
				{
					ix2 = 2*nmom - ix + 1 ;
				}
				if(iy<=nmom)
				{
					iy2 = iy - 1 ;
				}
				else
				{
					iy2 = 2*nmom - iy + 1 ;
				}
				Vq_temp[(2*ix-1)+(iy-1)*2*2*nmom]= Vq[iwb][ix2][iy2] ; 
				Vq_temp[2*ix+(iy-1)*2*2*nmom] = 0 ;
			}
		}
		fourn(Vq_temp,nmom2,2,-1) ;
		for(ix=0;ix<=nmom;ix++)
		{
			for(iy=0;iy<=nmom;iy++)	
			{

				ix2 = ix + 1 ;
				iy2 = iy + 1 ;
				Vq[iwb][ix][iy] = Vq_temp[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;	
			}
		}
	}
// fourier vq end


//	#pragma omp parallel for private(je)
	for(je=1;je<=Nener;je++) 
	{
		ener = enstep*(je-0.5) ;
		for(ix=0;ix<=nmom;ix++)
		{
			kx = step*(ix) ;
			for(iy=0;iy<=nmom;iy++)	
			{
				ky = step*(iy) ;

//				(n1+i*n2)/(de1-i*de2)

				ws1[je] = ener - (self_re[je][ix][iy] + self_imp_re[je][ix][iy]);
				xx1[je] = xi_band(kx,ky) + (renor_re[je][ix][iy] + renor_imp_re[je][ix][iy]) ;

				de_r[je] = pow(ws1[je],2) - pow(self_im[je][ix][iy]+self_imp_im[je][ix][iy],2) - pow(xx1[je],2) + pow(renor_im[je][ix][iy]+renor_imp_im[je][ix][iy],2) - pow(pairing_re[je][ix][iy]+pairing_imp_re[je][ix][iy],2) + pow(pairing_im[je][ix][iy]+pairing_imp_im[je][ix][iy],2) ;
				de_i[je] = 2*ws1[je]*(self_im[je][ix][iy]+self_imp_im[je][ix][iy]) + 2*xx1[je]*(renor_im[je][ix][iy]+renor_imp_im[je][ix][iy]) + 2*(pairing_re[je][ix][iy]+pairing_imp_re[je][ix][iy])*(pairing_im[je][ix][iy]+pairing_imp_im[je][ix][iy]) ;

				Gs_re[je][ix][iy] = (ws1[je] * de_r[je] + (self_im[je][ix][iy]+self_imp_im[je][ix][iy])*de_i[je])/(pow(de_r[je],2)+pow(de_i[je],2));
				Gx_re[je][ix][iy] = (xx1[je] * de_r[je] - (renor_im[je][ix][iy]+renor_imp_im[je][ix][iy])*de_i[je] )/(ow(de_r[je],2)+pow(de_i[je],2));
				Gp_re[je][ix][iy] = ((pairing_re[je][ix][iy]+pairing_imp_re[je][ix][iy]) * de_r[je] - (pairing_im[je][ix][iy]+pairing_imp_im[je][ix][iy])*de_i[je] )/(pow(de_r[je],2)+pow(de_i[je],2));

				Gs_im[je][ix][iy] = (ws1[je] * de_i[je] - (self_im[je][ix][iy]+self_imp_im[je][ix][iy])*de_r[je])/(pow(de_r[je],2)+pow(de_i[je],2));
				Gx_im[je][ix][iy] = (xx1[je] * de_i[je] + (renor_im[je][ix][iy]+renor_imp_im[je][ix][iy])*de_r[je] )/(pow(de_r[je],2)+pow(de_i[je],2));
				Gp_im[je][ix][iy] = ((pairing_re[je][ix][iy]+pairing_imp_re[je][ix][iy]) * de_i[je] + (pairing_im[je][ix][iy]+pairing_imp_im[je][ix][iy])*de_r[je] )/(pow(de_r[je],2)+pow(de_i[je],2));

				As[je][ix][iy] = - (1.0/M_PI) * Gs_im[je][ix][iy];
				Ax[je][ix][iy] = - (1.0/M_PI) * Gx_im[je][ix][iy];
				Ap[je][ix][iy] = - (1.0/M_PI) * Gp_im[je][ix][iy];

			}
		}
	}



// time
	time(&time_now);
	printf(ctime(&time_now));
	fprintf(fp,ctime(&time_now));
	running_time = time_now - time_0   ;
	printf("h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	
	fprintf(fp,"h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	


	error = 1 ;

	for(iter=1;iter<=N_iter;iter++)
	{
		printf("No of iterations is = %d .\t The error is = %.8f .\n" , iter,error);
		fprintf(fp,"No of iterations is = %d .\t The error is = %.8f .\n" , iter,error);

// fourier den start
		for(ie=1;ie<=Nener;ie++) 
		{
			for(ix=1;ix<=2*nmom;ix++)
			{
				for(iy=1;iy<=2*nmom;iy++)	
				{
					if(ix<=nmom)
					{
						ix2 = ix - 1 ;
					}
					else
					{
						ix2 = 2*nmom - ix + 1 ;
					}
					if(iy<=nmom)
					{
						iy2 = iy - 1 ;
					}
					else
					{
						iy2 = 2*nmom - iy + 1 ;
					}
					As_temp[(2*ix-1)+(iy-1)*2*2*nmom] = As[ie][ix2][iy2] ;
					As_temp[(2*ix)+(iy-1)*2*2*nmom] = 0 ;
					Ax_temp[(2*ix-1)+(iy-1)*2*2*nmom] = Ax[ie][ix2][iy2] ;
					Ax_temp[(2*ix)+(iy-1)*2*2*nmom] = 0 ;
					Ap_temp[(2*ix-1)+(iy-1)*2*2*nmom] = Ap[ie][ix2][iy2] ;
					Ap_temp[(2*ix)+(iy-1)*2*2*nmom] = 0 ;
				}
			}
			fourn(As_temp,nmom2,2,-1) ;
			fourn(Ax_temp,nmom2,2,-1) ;
			fourn(Ap_temp,nmom2,2,-1) ;
			for(ix=0;ix<=nmom;ix++)
			{
				for(iy=0;iy<=nmom;iy++)	
				{
					ix2 = ix + 1 ;
					iy2 = iy + 1 ;

					As[ie][ix][iy] = As_temp[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
					Ax[ie][ix][iy] = Ax_temp[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
					Ap[ie][ix][iy] = Ap_temp[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
				}
			}
		}
// fourier den end

		for(ix=0;ix<=nmom;ix++)
		{
			kx = step*(ix) ;

			for(iy=0;iy<=ix;iy++)	
			{
				ky = step*(iy) ;
				#pragma omp parallel for private(ie)
				for(ie=1;ie<=Nener;ie++) 
				{
					w = enstep*(ie-0.5) ;
					s_new_im[ie] = 0 ;
					x_new_im[ie] = 0 ;
					p_new_im[ie] = 0 ;

					for(iwb=iwbs;iwb<=iwbc;iwb++)
					{
						wboson = enstep*(iwb) ;
						if((ie-iwb>=0)&&(ie+iwb<=Nener))
						{
							s_new_im[ie] = s_new_im[ie]+enstep*Vq[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*As[ie-iwb][ix][iy] 
												-(fermi_d(w+wboson)+bose_d(+wboson))*As[ie+iwb][ix][iy]) ;	
							x_new_im[ie] = x_new_im[ie]+enstep*Vq[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*Ax[ie-iwb][ix][iy] 
												-(fermi_d(w+wboson)+bose_d(+wboson))*Ax[ie+iwb][ix][iy]) ;	
							p_new_im[ie] = p_new_im[ie]+enstep*Vq[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*Ap[ie-iwb][ix][iy] 
												-(fermi_d(w+wboson)+bose_d(+wboson))*Ap[ie+iwb][ix][iy]) ;	
						}
						else if((ie-iwb<0)&&(ie+iwb<=Nener))
						{
							s_new_im[ie] = s_new_im[ie]+enstep*Vq[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*As[iwb-ie][ix][iy] 
												-(fermi_d(w+wboson)+bose_d(+wboson))*As[ie+iwb][ix][iy]) ;	
							x_new_im[ie] = x_new_im[ie]+enstep*Vq[iwb][ix][iy]*(-(fermi_d(w-wboson)+bose_d(-wboson))*Ax[iwb-ie][ix][iy] 
												-(fermi_d(w+wboson)+bose_d(+wboson))*Ax[ie+iwb][ix][iy]) ;	
							p_new_im[ie] = p_new_im[ie]+enstep*Vq[iwb][ix][iy]*(-(fermi_d(w-wboson)+bose_d(-wboson))*Ap[iwb-ie][ix][iy] 
												-(fermi_d(w+wboson)+bose_d(+wboson))*Ap[ie+iwb][ix][iy]) ;	
						}
						else if((ie-iwb>=0)&&(ie+iwb>Nener))
						{
							s_new_im[ie] = s_new_im[ie]+enstep*Vq[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*As[ie-iwb][ix][iy] 
												-(fermi_d(w+wboson)+bose_d(+wboson))*As[Nener][ix][iy]) ;	
							x_new_im[ie] = x_new_im[ie]+enstep*Vq[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*Ax[ie-iwb][ix][iy] 
												-(fermi_d(w+wboson)+bose_d(+wboson))*Ax[Nener][ix][iy]) ;	
							p_new_im[ie] = p_new_im[ie]+enstep*Vq[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*Ap[ie-iwb][ix][iy] 
												-(fermi_d(w+wboson)+bose_d(+wboson))*Ap[Nener][ix][iy]) ;	
						}
					}

					self_im_new[ie][ix][iy] = M_PI*s_new_im[ie] ;
					self_im_new[ie][iy][ix] = M_PI*s_new_im[ie] ;
					renor_im_new[ie][ix][iy] = M_PI*x_new_im[ie];
					renor_im_new[ie][iy][ix] = M_PI*x_new_im[ie] ;
					pairing_im_new[ie][ix][iy] = M_PI*p_new_im[ie] ;
					pairing_im_new[ie][iy][ix] = -M_PI*p_new_im[ie] ;
//					pairing_im_new[ie][ix][iy] = 0 ;
//					pairing_im_new[ie][iy][ix] = 0 ;
				}

				#pragma omp parallel for private(ie)
				for(ie=1;ie<=Nener;ie++) 
				{
					s_new_re[ie] = 0 ;
					x_new_re[ie] = 0 ;
					p_new_re[ie] = 0 ;
					w = enstep*(ie-0.5) ;
					for(je=1;je<=Nener;je++) 
					{
						ener = enstep*(je-0.5)  ;
						aa = (pow(ener,2)-pow(w,2))/(pow(pow(ener,2)-pow(w,2),2)+pow(deltkk,2)) ;
						s_new_re[ie] = s_new_re[ie] + enstep*(1.0/M_PI) * (2*w*self_im_new[je][ix][iy]-2*w*self_im_new[ie][ix][iy]) * aa;
						x_new_re[ie] = x_new_re[ie] + enstep*(1.0/M_PI) * (2*ener*renor_im_new[je][ix][iy]-2*w*renor_im_new[ie][ix][iy]) * aa;
						p_new_re[ie] = p_new_re[ie] + enstep*(1.0/M_PI) * (2*ener*pairing_im_new[je][ix][iy]-2*w*pairing_im_new[ie][ix][iy]) * aa;
					}
					
					self_re_new[ie][ix][iy] = s_new_re[ie] ;
					self_re_new[ie][iy][ix] = s_new_re[ie] ;
					renor_re_new[ie][ix][iy] = x_new_re[ie] ;
					renor_re_new[ie][iy][ix] = x_new_re[ie] ;
					pairing_re_new[ie][ix][iy] = p_new_re[ie] ;
					pairing_re_new[ie][iy][ix] = -p_new_re[ie] ;
//					pairing_re_new[ie][ix][iy] = 0 ;
//					pairing_re_new[ie][iy][ix] = 0 ;
				}
			}
		}	


// inverse fourier
		for(ie=1;ie<=Nener;ie++) 
		{
			for(ix=1;ix<=2*nmom;ix++)
			{
				for(iy=1;iy<=2*nmom;iy++)
				{
					if(ix<=nmom)
					{
						ix2 = ix - 1 ;
					}
					else
					{
						ix2 = 2*nmom - ix + 1 ;
					}
					if(iy<=nmom)
					{
						iy2 = iy - 1 ;
					}
					else
					{
						iy2 = 2*nmom - iy + 1 ;
					}
					self_temp[(2*ix-1)+(iy-1)*2*2*nmom] = self_re_new[ie][ix2][iy2] ;
					self_temp[(2*ix)+(iy-1)*2*2*nmom] = self_im_new[ie][ix2][iy2] ;
					renor_temp[(2*ix-1)+(iy-1)*2*2*nmom] = renor_re_new[ie][ix2][iy2] ;
					renor_temp[(2*ix)+(iy-1)*2*2*nmom] = renor_im_new[ie][ix2][iy2] ;
					pairing_temp[(2*ix-1)+(iy-1)*2*2*nmom] = pairing_re_new[ie][ix2][iy2] ;
					pairing_temp[(2*ix)+(iy-1)*2*2*nmom] = pairing_im_new[ie][ix2][iy2] ;
				}
			}
			fourn(self_temp,nmom2,2,1) ;
			fourn(renor_temp,nmom2,2,1) ;
			fourn(pairing_temp,nmom2,2,1) ;
			for(ix=0;ix<=nmom;ix++)
			{
				for(iy=0;iy<=nmom;iy++)	
				{
					ix2 = ix + 1 ;
					iy2 = iy + 1 ;

					self_re_new[ie][ix][iy] = self_temp[(2*ix2-1)+(iy2-1)*2*2*nmom]  ;
					self_im_new[ie][ix][iy] = self_temp[(2*ix2)+(iy2-1)*2*2*nmom] ;
					renor_re_new[ie][ix][iy] = renor_temp[(2*ix2-1)+(iy2-1)*2*2*nmom]  ;
					renor_im_new[ie][ix][iy] = renor_temp[(2*ix2)+(iy2-1)*2*2*nmom] ;
					pairing_re_new[ie][ix][iy] = pairing_temp[(2*ix2-1)+(iy2-1)*2*2*nmom]  ;
					pairing_im_new[ie][ix][iy] = pairing_temp[(2*ix2)+(iy2-1)*2*2*nmom] ;
				}
			}
		}

		err_s = 0 ;
		err_x = 0 ;
		err_p = 0 ;
		
		for(ix=0;ix<=nmom;ix++)
		{
			for(iy=0;iy<=nmom;iy++)	
			{
				for(ie=1;ie<=Nener;ie++) 
				{
					err_s = err_s + pow(self_re_new[ie][ix][iy]-self_re[ie][ix][iy],2) + pow(self_im_new[ie][ix][iy]-self_im[ie][ix][iy],2) ;
					err_x = err_x + pow(renor_re_new[ie][ix][iy]-renor_re[ie][ix][iy],2) + pow(renor_im_new[ie][ix][iy]-renor_im[ie][ix][iy],2) ;
					err_p = err_p + pow(pairing_re_new[ie][ix][iy]-pairing_re[ie][ix][iy],2) + pow(pairing_im_new[ie][ix][iy]-pairing_im[ie][ix][iy],2) ;

				}
			}
		}
		
		error_s = sqrt(err_s) * pow(step,2) * enstep ;
		error_x = sqrt(err_x) * pow(step,2) * enstep ;
		error_p = sqrt(err_p) * pow(step,2) * enstep ;
		error = error_s + error_x + error_p ;

		printf("error : %.8f\t%.8f\t%.8f\n", error_s,error_x,error_p) ;

		mixing = 0.9 ;
		for(ix=0;ix<=nmom;ix++)
		{
			for(iy=0;iy<=nmom;iy++)	
			{
				for(ie=1;ie<=Nener;ie++) 
				{

					self_re[ie][ix][iy] = (1.0-mixing)*self_re[ie][ix][iy] + mixing*self_re_new[ie][ix][iy] ;
					self_im[ie][ix][iy] = (1.0-mixing)*self_im[ie][ix][iy] + mixing*self_im_new[ie][ix][iy] ;
					renor_re[ie][ix][iy] = (1.0-mixing)*renor_re[ie][ix][iy] + mixing*renor_re_new[ie][ix][iy] ;
					renor_im[ie][ix][iy] = (1.0-mixing)*renor_im[ie][ix][iy] + mixing*renor_im_new[ie][ix][iy] ;
					pairing_re[ie][ix][iy] = (1.0-mixing)*pairing_re[ie][ix][iy] + mixing*pairing_re_new[ie][ix][iy] ;
					pairing_im[ie][ix][iy] = (1.0-mixing)*pairing_im[ie][ix][iy] + mixing*pairing_im_new[ie][ix][iy] ;

				}
			}
		}		

	// time
		time(&time_now);
		printf(ctime(&time_now));
		fprintf(fp,ctime(&time_now));
		running_time = time_now - time_0   ;
		printf("h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	
		fprintf(fp,"h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	



		#pragma omp parallel for private(je)
		for(je=1;je<=Nener;je++) 
		{
			ener = enstep*(je-0.5) ;
			for(ix=0;ix<=nmom;ix++)
			{
				kx = step*(ix) ;
				for(iy=0;iy<=nmom;iy++)	
				{
					ky = step*(iy) ;

//					(n1+i*n2)/(de1-i*de2)

					ws1[je] = ener - (self_re[je][ix][iy] + self_imp_re[je][ix][iy]);
					xx1[je] = xi_band(kx,ky) + (renor_re[je][ix][iy] + renor_imp_re[je][ix][iy]) ;

					de_r[je] = pow(ws1[je],2) - pow(self_im[je][ix][iy]+self_imp_im[je][ix][iy],2) - pow(xx1[je],2) + pow(renor_im[je][ix][iy]+renor_imp_im[je][ix][iy],2) - pow(pairing_re[je][ix][iy]+pairing_imp_re[je][ix][iy],2) + pow(pairing_im[je][ix][iy]+pairing_imp_im[je][ix][iy],2) ;
					de_i[je] = 2*ws1[je]*(self_im[je][ix][iy]+self_imp_im[je][ix][iy]) + 2*xx1[je]*(renor_im[je][ix][iy]+renor_imp_im[je][ix][iy]) + 2*(pairing_re[je][ix][iy]+pairing_imp_re[je][ix][iy])*(pairing_im[je][ix][iy]+pairing_imp_im[je][ix][iy]) ;

					Gs_re[je][ix][iy] = (ws1[je] * de_r[je] + (self_im[je][ix][iy]+self_imp_im[je][ix][iy])*de_i[je])/(pow(de_r[je],2)+pow(de_i[je],2));
					Gx_re[je][ix][iy] = (xx1[je] * de_r[je] - (renor_im[je][ix][iy]+renor_imp_im[je][ix][iy])*de_i[je] )/(pow(de_r[je],2)+pow(de_i[je],2));
					Gp_re[je][ix][iy] = ((pairing_re[je][ix][iy]+pairing_imp_re[je][ix][iy]) * de_r[je] - (pairing_im[je][ix][iy]+pairing_imp_im[je][ix][iy])*de_i[je] )/(pow(de_r[je],2)+pow(de_i[je],2));

					Gs_im[je][ix][iy] = (ws1[je] * de_i[je] - (self_im[je][ix][iy]+self_imp_im[je][ix][iy])*de_r[je])/(pow(de_r[je],2)+pow(de_i[je],2));
					Gx_im[je][ix][iy] = (xx1[je] * de_i[je] + (renor_im[je][ix][iy]+renor_imp_im[je][ix][iy])*de_r[je] )/(pow(de_r[je],2)+pow(de_i[je],2));
					Gp_im[je][ix][iy] = ((pairing_re[je][ix][iy]+pairing_imp_re[je][ix][iy]) * de_i[je] + (pairing_im[je][ix][iy]+pairing_imp_im[je][ix][iy])*de_r[je] )/(pow(de_r[je],2)+pow(de_i[je],2));

					As[je][ix][iy] = - (1.0/M_PI) * Gs_im[je][ix][iy];
					Ax[je][ix][iy] = - (1.0/M_PI) * Gx_im[je][ix][iy];
					Ap[je][ix][iy] = - (1.0/M_PI) * Gp_im[je][ix][iy];

				}
			}
		}



		printf("iteration 3\n");
		fprintf(fp,"iteration 3\n" );	

		error3 = 1 ;

		for(iter3=1;iter3<=N_iter3;iter3++)
		{
			printf("No of iterations is = %d .\t The error is = %.8f .\n" , iter3,error3);
			fprintf(fp,"No of iterations is = %d .\t The error is = %.8f .\n" , iter3,error3);

	// fourier den start
			for(ie=1;ie<=Nener;ie++) 
			{
				for(ix=1;ix<=2*nmom;ix++)
				{
					for(iy=1;iy<=2*nmom;iy++)
					{
						if(ix<=nmom)
						{
							ix2 = ix - 1 ;
						}
						else
						{
							ix2 = 2*nmom - ix + 1 ;
						}
						if(iy<=nmom)
						{
							iy2 = iy - 1 ;
						}
						else
						{
							iy2 = 2*nmom - iy + 1 ;
						}
						Gs_temp[(2*ix-1)+(iy-1)*2*2*nmom] = Gs_re[ie][ix2][iy2] ;
						Gs_temp[(2*ix)+(iy-1)*2*2*nmom] = Gs_im[ie][ix2][iy2] ;
						Gx_temp[(2*ix-1)+(iy-1)*2*2*nmom] = Gx_re[ie][ix2][iy2] ;
						Gx_temp[(2*ix)+(iy-1)*2*2*nmom] = Gx_im[ie][ix2][iy2] ;
						Gp_temp[(2*ix-1)+(iy-1)*2*2*nmom] = Gp_re[ie][ix2][iy2] ;
						Gp_temp[(2*ix)+(iy-1)*2*2*nmom] = Gp_im[ie][ix2][iy2] ;
					}
				}
				fourn(Gs_temp,nmom2,2,-1) ;
				fourn(Gx_temp,nmom2,2,-1) ;
				fourn(Gp_temp,nmom2,2,-1) ;
				for(ix=0;ix<=nmom;ix++)
				{
					for(iy=0;iy<=nmom;iy++)	
					{
						ix2 = ix + 1 ;
						iy2 = iy + 1 ;

						Gs_re[ie][ix][iy] = Gs_temp[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
						Gs_im[ie][ix][iy] = Gs_temp[(2*ix2)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
						Gx_re[ie][ix][iy] = Gx_temp[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
						Gx_im[ie][ix][iy] = Gx_temp[(2*ix2)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
						Gp_re[ie][ix][iy] = Gp_temp[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
						Gp_im[ie][ix][iy] = Gp_temp[(2*ix2)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
					}
				}
			}
	// fourier den end

			printf("Fourier\n");
			fprintf(fp,"Fourier\n");

		// time
			time(&time_now);
			printf(ctime(&time_now));
			fprintf(fp,ctime(&time_now));
			running_time = time_now - time_0   ;
			printf("h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	
			fprintf(fp,"h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	


			for(ix=0;ix<=nmom;ix++)
			{
				kx = step*(ix) ;
				for(iy=0;iy<=nmom;iy++)	
				{
					ky = step*(iy) ;
					#pragma omp parallel for private(ie)
					for(ie=1;ie<=Nener;ie++) 
					{
						w = enstep*(ie-0.5) ;
						self_imp_re_new[ie][ix][iy] = Vq_imp[ix][iy]*Gs_re[ie][ix][iy] ;
						self_imp_im_new[ie][ix][iy] = Vq_imp[ix][iy]*Gs_im[ie][ix][iy] ;
						renor_imp_re_new[ie][ix][iy] = Vq_imp[ix][iy]*Gx_re[ie][ix][iy] ;
						renor_imp_im_new[ie][ix][iy] = Vq_imp[ix][iy]*Gx_im[ie][ix][iy] ;
						pairing_imp_re_new[ie][ix][iy] = -Vq_imp[ix][iy]*Gp_re[ie][ix][iy] ;
						pairing_imp_im_new[ie][ix][iy] = -Vq_imp[ix][iy]*Gp_im[ie][ix][iy] ;

					}
				}
			}

			printf("imp\n");
			fprintf(fp,"imp\n");
		// time
			time(&time_now);
			printf(ctime(&time_now));
			fprintf(fp,ctime(&time_now));
			running_time = time_now - time_0   ;
			printf("h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	
			fprintf(fp,"h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	


	// inverse fourier
			for(ie=1;ie<=Nener;ie++) 
			{
				for(ix=1;ix<=2*nmom;ix++)
				{
					for(iy=1;iy<=2*nmom;iy++)
					{
						if(ix<=nmom)
						{
							ix2 = ix - 1 ;
						}
						else
						{
							ix2 = 2*nmom - ix + 1 ;
						}
						if(iy<=nmom)
						{
							iy2 = iy - 1 ;
						}
						else
						{
							iy2 = 2*nmom - iy + 1 ;
						}
						self_temp[(2*ix-1)+(iy-1)*2*2*nmom] = self_imp_re_new[ie][ix2][iy2] ;
						self_temp[(2*ix)+(iy-1)*2*2*nmom] = self_imp_im_new[ie][ix2][iy2] ;
						renor_temp[(2*ix-1)+(iy-1)*2*2*nmom] = renor_imp_re_new[ie][ix2][iy2] ;
						renor_temp[(2*ix)+(iy-1)*2*2*nmom] = renor_imp_im_new[ie][ix2][iy2] ;
						pairing_temp[(2*ix-1)+(iy-1)*2*2*nmom] = pairing_imp_re_new[ie][ix2][iy2] ;
						pairing_temp[(2*ix)+(iy-1)*2*2*nmom] = pairing_imp_im_new[ie][ix2][iy2] ;
					}
				}
				fourn(self_temp,nmom2,2,1) ;
				fourn(renor_temp,nmom2,2,1) ;
				fourn(pairing_temp,nmom2,2,1) ;
				for(ix=0;ix<=nmom;ix++)
				{
					for(iy=0;iy<=nmom;iy++)	
					{
						ix2 = ix + 1 ;
						iy2 = iy + 1 ;

						self_imp_re_new[ie][ix][iy] = self_temp[(2*ix2-1)+(iy2-1)*2*2*nmom]  ;
						self_imp_im_new[ie][ix][iy] = self_temp[(2*ix2)+(iy2-1)*2*2*nmom] ;
						renor_imp_re_new[ie][ix][iy] = renor_temp[(2*ix2-1)+(iy2-1)*2*2*nmom]  ;
						renor_imp_im_new[ie][ix][iy] = renor_temp[(2*ix2)+(iy2-1)*2*2*nmom] ;
						pairing_imp_re_new[ie][ix][iy] = pairing_temp[(2*ix2-1)+(iy2-1)*2*2*nmom]  ;
						pairing_imp_im_new[ie][ix][iy] = pairing_temp[(2*ix2)+(iy2-1)*2*2*nmom] ;
					}
				}
			}

			printf("Fourier\n");
			fprintf(fp,"Fourier\n");
		// time
			time(&time_now);
			printf(ctime(&time_now));
			fprintf(fp,ctime(&time_now));
			running_time = time_now - time_0   ;
			printf("h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	
			fprintf(fp,"h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	

			err_s = 0 ;
			err_x = 0 ;
			err_p = 0 ;
		
			for(ix=0;ix<=nmom;ix++)
			{
				for(iy=0;iy<=nmom;iy++)	
				{
					for(ie=1;ie<=Nener;ie++) 
					{
						err_s = err_s + pow(self_imp_re_new[ie][ix][iy]-self_imp_re[ie][ix][iy],2) + pow(self_imp_im_new[ie][ix][iy]-self_imp_im[ie][ix][iy],2) ;
						err_x = err_x + pow(renor_imp_re_new[ie][ix][iy]-renor_imp_re[ie][ix][iy],2) + pow(renor_imp_im_new[ie][ix][iy]-renor_imp_im[ie][ix][iy],2) ;
						err_p = err_p + pow(pairing_imp_re_new[ie][ix][iy]-pairing_imp_re[ie][ix][iy],2) + pow(pairing_imp_im_new[ie][ix][iy]-pairing_imp_im[ie][ix][iy],2) ;
					}
				}
			}
		
			error_s = sqrt(err_s) * pow(step,2) * enstep ;
			error_x = sqrt(err_x) * pow(step,2) * enstep ;
			error_p = sqrt(err_p) * pow(step,2) * enstep ;
			error3 = error_s + error_x + error_p ;

			printf("error : %.8f\t%.8f\t%.8f\n", error_s,error_x,error_p) ;
			fprintf(fp,"error : %.8f\t%.8f\t%.8f\n", error_s,error_x,error_p) ;

			mixing3 = 1.0 ;
			for(ix=0;ix<=nmom;ix++)
			{
				for(iy=0;iy<=nmom;iy++)	
				{
					#pragma omp parallel for private(ie)
					for(ie=1;ie<=Nener;ie++) 
					{
						self_imp_re[ie][ix][iy] = (1.0-mixing3)*self_imp_re[ie][ix][iy] + mixing3*self_imp_re_new[ie][ix][iy] ;
						self_imp_im[ie][ix][iy] = (1.0-mixing3)*self_imp_im[ie][ix][iy] + mixing3*self_imp_im_new[ie][ix][iy] ;
						renor_imp_re[ie][ix][iy] = (1.0-mixing3)*renor_imp_re[ie][ix][iy] + mixing3*renor_imp_re_new[ie][ix][iy] ;
						renor_imp_im[ie][ix][iy] = (1.0-mixing3)*renor_imp_im[ie][ix][iy] + mixing3*renor_imp_im_new[ie][ix][iy] ;
						pairing_imp_re[ie][ix][iy] = (1.0-mixing3)*pairing_imp_re[ie][ix][iy] + mixing3*pairing_imp_re_new[ie][ix][iy] ;
						pairing_imp_im[ie][ix][iy] = (1.0-mixing3)*pairing_imp_im[ie][ix][iy] + mixing3*pairing_imp_im_new[ie][ix][iy] ;
					}
				}
			}

	//		#pragma omp parallel for private(je)
			for(je=1;je<=Nener;je++) 
			{
				ener = enstep*(je) ;
				for(ix=0;ix<=nmom;ix++)
				{
					kx = step*(ix) ;
					for(iy=0;iy<=nmom;iy++)	
					{
						ky = step*(iy) ;

		//				(n1+i*n2)/(de1-i*de2)

						ws1[je] = ener - (self_re[je][ix][iy] + self_imp_re[je][ix][iy]);
						xx1[je] = xi_band(kx,ky) + (renor_re[je][ix][iy] + renor_imp_re[je][ix][iy]) ;

						de_r[je] = pow(ws1[je],2) - pow(self_im[je][ix][iy]+self_imp_im[je][ix][iy],2) - pow(xx1[je],2) + pow(renor_im[je][ix][iy]+renor_imp_im[je][ix][iy],2) - pow(pairing_re[je][ix][iy]+pairing_imp_re[je][ix][iy],2) + pow(pairing_im[je][ix][iy]+pairing_imp_im[je][ix][iy],2) ;
						de_i[je] = 2*ws1[je]*(self_im[je][ix][iy]+self_imp_im[je][ix][iy]) + 2*xx1[je]*(renor_im[je][ix][iy]+renor_imp_im[je][ix][iy]) + 2*(pairing_re[je][ix][iy]+pairing_imp_re[je][ix][iy])*(pairing_im[je][ix][iy]+pairing_imp_im[je][ix][iy]) ;

						Gs_re[je][ix][iy] = (ws1[je] * de_r[je] + (self_im[je][ix][iy]+self_imp_im[je][ix][iy])*de_i[je])/(pow(de_r[je],2)+pow(de_i[je],2));
						Gx_re[je][ix][iy] = (xx1[je] * de_r[je] - (renor_im[je][ix][iy]+renor_imp_im[je][ix][iy])*de_i[je] )/(pow(de_r[je],2)+pow(de_i[je],2));
						Gp_re[je][ix][iy] = ((pairing_re[je][ix][iy]+pairing_imp_re[je][ix][iy]) * de_r[je] - (pairing_im[je][ix][iy]+pairing_imp_im[je][ix][iy])*de_i[je] )/(pow(de_r[je],2)+pow(de_i[je],2));
				
						Gs_im[je][ix][iy] = (ws1[je] * de_i[je] - (self_im[je][ix][iy]+self_imp_im[je][ix][iy])*de_r[je])/(pow(de_r[je],2)+pow(de_i[je],2));
						Gx_im[je][ix][iy] = (xx1[je] * de_i[je] + (renor_im[je][ix][iy]+renor_imp_im[je][ix][iy])*de_r[je] )/(pow(de_r[je],2)+pow(de_i[je],2));
						Gp_im[je][ix][iy] = ((pairing_re[je][ix][iy]+pairing_imp_re[je][ix][iy]) * de_i[je] + (pairing_im[je][ix][iy]+pairing_imp_im[je][ix][iy])*de_r[je] )/(pow(de_r[je],2)+pow(de_i[je],2));
				
						As[je][ix][iy] = - (1.0/M_PI) * Gs_im[je][ix][iy];
						Ax[je][ix][iy] = - (1.0/M_PI) * Gx_im[je][ix][iy];
						Ap[je][ix][iy] = - (1.0/M_PI) * Gp_im[je][ix][iy];

					}
				}
			}

		// time
			time(&time_now);
			printf(ctime(&time_now));
			fprintf(fp,ctime(&time_now));
			running_time = time_now - time_0   ;
			printf("h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	
			fprintf(fp,"h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	

			if(error3<=tol3)
			{
				iter3 = N_iter3 + 101 ;
			}
		
		}



	// time
		time(&time_now);
		printf(ctime(&time_now));
		fprintf(fp,ctime(&time_now));
		running_time = time_now - time_0   ;
		printf("h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	
		fprintf(fp,"h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	

		if(error<=tol)
		{
			iter = N_iter + 101 ;
		}
		
	}





//       ppppp
	#pragma omp parallel for private(je)
	for(je=1;je<=Nener;je++) 
	{
		ener = enstep*(je-0.5) ;
		Ns_re[je] = 0 ;
		Nx_re[je] = 0 ;
		
		for(ix=0;ix<=nmom;ix++)
		{
			kx = step*(ix) ;
			for(iy=0;iy<=nmom;iy++)	
			{
				ky = step*(iy) ;
				if((ix==0)&&(iy==0))
				{
					Ns_re[je] = Ns_re[je] + step*step*As[je][ix][iy]/2.0 ;
					Nx_re[je] = Nx_re[je] + step*step*Ax[je][ix][iy]/2.0 ;
					Ns_im[je] = Ns_im[je] - M_PI*step*step*Gs_re[je][ix][iy]/2.0 ;
					Nx_im[je] = Nx_im[je] - M_PI*step*step*Gx_re[je][ix][iy]/2.0 ;
				}
				else if((ix==nmom)&&(iy==nmom))
				{
					Ns_re[je] = Ns_re[je] + step*step*As[je][ix][iy]/2.0 ;
					Nx_re[je] = Nx_re[je] + step*step*Ax[je][ix][iy]/2.0 ;
					Ns_im[je] = Ns_im[je] - M_PI*step*step*Gs_re[je][ix][iy]/2.0 ;
					Nx_im[je] = Nx_im[je] - M_PI*step*step*Gx_re[je][ix][iy]/2.0 ;
				}
				else if((ix==0)&&(iy==nmom))
				{
					Ns_re[je] = Ns_re[je] + step*step*As[je][ix][iy]/2.0 ;
					Nx_re[je] = Nx_re[je] + step*step*Ax[je][ix][iy]/2.0 ;
					Ns_im[je] = Ns_im[je] - M_PI*step*step*Gs_re[je][ix][iy]/2.0 ;
					Nx_im[je] = Nx_im[je] - M_PI*step*step*Gx_re[je][ix][iy]/2.0 ;
				}
				else if((ix==nmom)&&(iy==0))
				{
					Ns_re[je] = Ns_re[je] + step*step*As[je][ix][iy]/2.0 ;
					Nx_re[je] = Nx_re[je] + step*step*Ax[je][ix][iy]/2.0 ;
					Ns_im[je] = Ns_im[je] - M_PI*step*step*Gs_re[je][ix][iy]/2.0 ;
					Nx_im[je] = Nx_im[je] - M_PI*step*step*Gx_re[je][ix][iy]/2.0 ;
				}
				else if((ix==0)&&(iy!=0))
				{
					Ns_re[je] = Ns_re[je] + 2*step*step*As[je][ix][iy]/2.0 ;
					Nx_re[je] = Nx_re[je] + 2*step*step*Ax[je][ix][iy]/2.0 ;
					Ns_im[je] = Ns_im[je] - 2*M_PI*step*step*Gs_re[je][ix][iy]/2.0 ;
					Nx_im[je] = Nx_im[je] - 2*M_PI*step*step*Gx_re[je][ix][iy]/2.0 ;				}
				else if((ix!=0)&&(iy==0))
				{
					Ns_re[je] = Ns_re[je] + 2*step*step*As[je][ix][iy]/2.0 ;
					Nx_re[je] = Nx_re[je] + 2*step*step*Ax[je][ix][iy]/2.0 ;
					Ns_im[je] = Ns_im[je] - 2*M_PI*step*step*Gs_re[je][ix][iy]/2.0 ;
					Nx_im[je] = Nx_im[je] - 2*M_PI*step*step*Gx_re[je][ix][iy]/2.0 ;
				}
				else if((ix==nmom)&&(iy!=nmom))
				{
					Ns_re[je] = Ns_re[je] + 2*step*step*As[je][ix][iy]/2.0 ;
					Nx_re[je] = Nx_re[je] + 2*step*step*Ax[je][ix][iy]/2.0 ;
					Ns_im[je] = Ns_im[je] - 2*M_PI*step*step*Gs_re[je][ix][iy]/2.0 ;
					Nx_im[je] = Nx_im[je] - 2*M_PI*step*step*Gx_re[je][ix][iy]/2.0 ;
				}
				else if((ix!=nmom)&&(iy==nmom))
				{
					Ns_re[je] = Ns_re[je] + 2*step*step*As[je][ix][iy]/2.0 ;
					Nx_re[je] = Nx_re[je] + 2*step*step*Ax[je][ix][iy]/2.0 ;
					Ns_im[je] = Ns_im[je] - 2*M_PI*step*step*Gs_re[je][ix][iy]/2.0 ;
					Nx_im[je] = Nx_im[je] - 2*M_PI*step*step*Gx_re[je][ix][iy]/2.0 ;
				}
				else
				{
					Ns_re[je] = Ns_re[je] + 4*step*step*As[je][ix][iy]/2.0 ;
					Nx_re[je] = Nx_re[je] + 4*step*step*Ax[je][ix][iy]/2.0 ;
					Ns_im[je] = Ns_im[je] - 2*M_PI*step*step*Gs_re[je][ix][iy]/2.0 ;
					Nx_im[je] = Nx_im[je] - 2*M_PI*step*step*Gx_re[je][ix][iy]/2.0 ;
				}
			}
		}
	}


	fp7 = fopen("d7.d","w");

	for(je=Nener;je>=1;je--) 
	{	
		ener = -enstep*(je-0.5) ;
		fprintf(fp7,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",ener,Ns_re[je],-Ns_im[je],-Nx_re[je],Nx_im[je],(Ns_re[je] - Nx_re[je]),(-Ns_im[je] + Nx_im[je])) ;
	}
	for(je=1;je<=Nener;je++) 
	{	
		ener = enstep*(je-0.5) ;
		fprintf(fp7,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",ener,Ns_re[je],Ns_im[je],Nx_re[je],Nx_im[je],(Ns_re[je] + Nx_re[je]),(Ns_im[je] + Nx_im[je])) ;
	}

	fclose(fp7) ;



	fpmi = fopen("Aminus.d","w") ;
	Ami = 0 ;
	for(je=1;je<=Nener;je++) 
	{
		ener = enstep*(je-0.5) ;
		for(ix=0;ix<=nmom;ix++)
		{
			kx = step*(ix) ;
			for(iy=0;iy<=nmom;iy++)	
			{
				ky = step*(iy) ;
				if(As[je][ix][iy]+Ax[je][ix][iy]<0)
				{
					fprintf(fpmi,"%f\t%f\t%f\t%f\n",kx,ky,ener,As[je][ix][iy]+Ax[je][ix][iy]) ;
					Ami = Ami + 1 ;
				}
				if(As[je][ix][iy]+Ax[je][ix][iy]<0)
				{
					fprintf(fpmi,"%f\t%f\t%f\t%f\n",kx,ky,-ener,As[je][ix][iy]-Ax[je][ix][iy]) ;
					Ami = Ami + 1 ;
				}
			}
		}
	}
	printf("number of the minus A : %d\n", Ami) ;
	fprintf(fp,"number of the minus A : %d\n", Ami) ;
	fclose(fpmi) ;




	fp2 = fopen("A000.dat","w");
	for(ix=0;ix<=nmom;ix++)	
	{
		kx = step*(ix) ;

		for(iy=0;iy<=nmom;iy++)	
		{
			ky = step*(iy) ;
//			if((ix%2==0)&&(iy%2==0))
//			{
				den0 = As[1][ix][iy] + Ax[1][ix][iy] ;
				fprintf(fp2,"%f\t%f\t%f\n" ,kx,ky,den0) ;
//			}
		}
	}
	fclose(fp2);


	sprintf(fid_mdc , "mdc");

	mkdir(fid_mdc,1017);

	for(ie=1;ie<=50;ie++)
	{
		w = (ie-0.5)*enstep ;
		sprintf(fi_mdc, "mdc/p%d.dat",ie) ;
		fp200 = fopen(fi_mdc,"w") ;
		for(ix=0;ix<=nmom;ix++)	
		{
			kx = step*(ix) ;

			for(iy=0;iy<=nmom;iy++)	
			{
				ky = step*(iy) ;
				den0 = As[ie][ix][iy] + Ax[ie][ix][iy] ;
				fprintf(fp200,"%f\t%f\t%f\n" ,kx,ky,den0) ;
			}
		}
		fclose(fp200);
	}
	for(ie=1;ie<=50;ie++)
	{
		w = -(ie-0.5)*enstep ;
		sprintf(fi_mdc, "mdc/m%d.dat",ie) ;
		fp200 = fopen(fi_mdc,"w") ;
		for(ix=0;ix<=nmom;ix++)	
		{
			kx = step*(ix) ;
			for(iy=0;iy<=nmom;iy++)	
			{
				ky = step*(iy) ;
				den0 = As[ie][ix][iy] - Ax[ie][ix][iy] ;
				fprintf(fp200,"%f\t%f\t%f\n" , kx,ky,den0) ;
			}
		}
		fclose(fp200);
	}

	fp6 = fopen("mat6.d","w");

	for (ix = 0; ix <= nmom; ix++)
	{
		kx = step*(ix) ;
		band_0[ix] = 0 ;
		min_band = 430 ;
		
		for(iy=0;iy<=nmom;iy++)	
		{
			ky = step*(iy) ;
		
			if(min_band>=fabs(xi_band(kx,ky)+ renor_re[0][ix][iy]))
			{
				min_band = fabs(xi_band(kx,ky)+ renor_re[0][ix][iy]) ;
				band_0[ix] = iy ;
			}
		}
		ky0 = step*(band_0[ix]) ;
		fprintf(fp6,"%d\t%d\t%f\t%f\t%f\n",ix,band_0[ix],kx,ky0,180*atan((1.0-ky0)/(1.0-kx))/M_PI);
	}
	fclose(fp6);


	
	sprintf(fid4 , "dat3");

	mkdir(fid4,1017);
	for(ix=0;ix<=nmom;ix++)
	{

		kx = step*(ix) ;
		iy = band_0[ix] ;
		ky = step*(iy) ;
		
		sprintf(fi_se4, "dat3/%d.dat",ix) ;
		fp31 = fopen(fi_se4,"w") ;
		for(ie=1;ie<=Nener;ie++)
		{
			w = enstep*(ie-0.5) ;
			fprintf(fp31, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",w,self_re[ie][ix][iy],self_im[ie][ix][iy],renor_re[ie][ix][iy],renor_im[ie][ix][iy],pairing_re[ie][ix][iy],pairing_im[ie][ix][iy],As[ie][ix][iy],Ax[ie][ix][iy],Ap[ie][ix][iy]) ;

		}

		fclose(fp31) ;
	}

	
	sprintf(fid5 , "dat4");

	mkdir(fid5,1017);
	for(ix=0;ix<=nmom;ix++)
	{

		kx = step*(ix) ;
		iy = 0 ;
		ky = step*(iy) ;
		
		sprintf(fi_se5, "dat4/%d.dat",ix) ;
		fp35 = fopen(fi_se5,"w") ;
		for(ie=0;ie<=Nener;ie++)
		{
			w = enstep*(ie) ;
			fprintf(fp35, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",w,self_re[ie][ix][iy],self_im[ie][ix][iy],renor_re[ie][ix][iy],renor_im[ie][ix][iy],pairing_re[ie][ix][iy],pairing_im[ie][ix][iy],As[ie][ix][iy],Ax[ie][ix][iy],Ap[ie][ix][iy]) ;
		}
		fclose(fp35) ;
	}


	sprintf(fid_angle , "sxp");
	sprintf(fid_edc , "edc3");
	mkdir(fid_angle,1017);
	mkdir(fid_edc,1017);
	for(itheta=0;itheta<=9;itheta++)
	{
		theta = itheta * 5;
		sprintf(fi_edc, "edc3/%d.dat",theta) ;

		fp_edc = fopen(fi_edc,"w") ;

		for(ix=0;ix<=nmom5;ix++)
		{
			kx = ix*step2 ;
			ky = 1.0-(1.0-kx)*tan(M_PI*theta/180) ;
			ix1 = 0 ;
			iy1 = 0 ;
			count_x = 0 ;
			count_y = 0 ;

			for(ix2=0;ix2<=nmom-1;ix2++)
			{
				kx1 = ix2*step;
				kx3 = (ix2+1)*step;
				if((kx1<=fabs(kx))&&(fabs(kx)<=kx3)&&(count_x==0))
				{
					ix1 = ix2;
					count_x = count_x + 1 ;
				}
			}

			if(count_x==0)
			{	
				printf("ix =%d\n",ix) ;
			}

			for(iy2=0;iy2<=nmom-1;iy2++)
			{
				ky1 = iy2*step;
				ky3 = (iy2+1)*step;
				if((ky1<=fabs(ky))&&(fabs(ky)<=ky3)&&(count_y==0))
				{
					iy1 = iy2;
					count_y = count_y + 1 ;
				}
			}
			if(count_y==0)
			{	
				printf("iy =%d\n",iy) ;
			}

			ix3 = ix1+1;
			iy3 = iy1+1;
			kx1 = ix1*step;
			kx3 = (ix1+1)*step;
			ky1 = iy1*step;
			ky3 = (iy1+1)*step;

			for(ie=Nener;ie>=1;ie--)
			{
				w = -enstep*(ie-0.5) ;
				iee = Nener-ie + 1 ; 
				A11 = As[ie][ix1][iy1] - Ax[ie][ix1][iy1] ;
				A31 = As[ie][ix3][iy1] - Ax[ie][ix3][iy1] ;
				A13 = As[ie][ix1][iy3] - Ax[ie][ix1][iy3] ;
				A33 = As[ie][ix3][iy3] - Ax[ie][ix3][iy3] ;
				Akxky[iee][ix] = ((kx3-fabs(kx))*(ky3-fabs(ky))* A11+(fabs(kx)-kx1)*(ky3-fabs(ky))* A31+(kx3-fabs(kx))*(fabs(ky)-ky1)* A13+(fabs(kx)-kx1)*(fabs(ky)-ky1)* A33)/((kx3-kx1)*(kx3-kx1));
				fprintf(fp_edc,"%f\t%f\t%f\t%f\n", w, kx, ky, Akxky[iee][ix]) ;
			}
			for(ie=1;ie<=Nener;ie++)
			{
				w = enstep*(ie-0.5) ;
				iee = Nener + ie + 1 ;
				A11 = As[ie][ix1][iy1] + Ax[ie][ix1][iy1] ;
				A31 = As[ie][ix3][iy1] + Ax[ie][ix3][iy1] ;
				A13 = As[ie][ix1][iy3] + Ax[ie][ix1][iy3] ;
				A33 = As[ie][ix3][iy3] + Ax[ie][ix3][iy3] ;
				Akxky[iee][ix] = ((kx3-fabs(kx))*(ky3-fabs(ky))* A11+(fabs(kx)-kx1)*(ky3-fabs(ky))* A31+(kx3-fabs(kx))*(fabs(ky)-ky1)* A13+(fabs(kx)-kx1)*(fabs(ky)-ky1)* A33)/((kx3-kx1)*(kx3-kx1));
				fprintf(fp_edc,"%f\t%f\t%f\t%f\n", w, kx, ky, Akxky[iee][ix]) ;
			}
		}
		fclose(fp_edc) ;
		
		sprintf(fi_angle, "sxp/%d.dat",theta) ;
		fp_angle = fopen(fi_angle,"w") ;
		for(ie=Nener;ie>=1;ie--)
		{
			w = -enstep*(ie-0.5) ;
			iee = Nener-ie + 1 ; 
			temp = 0;
			for(ix=0;ix<=nmom5;ix++)
			{
				if(temp<Akxky[iee][ix])
				{
					temp = Akxky[iee][ix];
					ix2 = ix ;
				}
			}
			kx = ix2*step2 ;
			ky = 1.0-(1.0-kx)*tan(M_PI*theta/180) ;
			ix1 = 0 ;
			iy1 = 0 ;
			count_x = 0 ;
			count_y = 0 ;

			for(ix2=0;ix2<=nmom-1;ix2++)
			{
				kx1 = ix2*step;
				kx3 = (ix2+1)*step;
				if((kx1<=fabs(kx))&&(fabs(kx)<=kx3)&&(count_x==0))
				{
					ix1 = ix2;
					count_x = count_x + 1 ;
				}
			}
			for(iy2=0;iy2<=nmom-1;iy2++)
			{
				ky1 = iy2*step;
				ky3 = (iy2+1)*step;
				if((ky1<=fabs(ky))&&(fabs(ky)<=ky3)&&(count_y==0))
				{
					iy1 = iy2;
					count_y = count_y + 1 ;
				}
			}
			ix3 = ix1+1;
			iy3 = iy1+1;
			kx1 = ix1*step;
			kx3 = (ix1+1)*step;
			ky1 = iy1*step;
			ky3 = (iy1+1)*step;

			s1kxky = -((kx3-fabs(kx))*(ky3-fabs(ky))* self_re[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* self_re[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* self_re[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* self_re[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			s2kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* self_im[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* self_im[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* self_im[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* self_im[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			x1kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* renor_re[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* renor_re[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* renor_re[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* renor_re[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			x2kxky = -((kx3-fabs(kx))*(ky3-fabs(ky))* renor_im[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* renor_im[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* renor_im[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* renor_im[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			p1kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* pairing_re[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* pairing_re[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* pairing_re[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* pairing_re[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			p2kxky = -((kx3-fabs(kx))*(ky3-fabs(ky))* pairing_im[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* pairing_im[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* pairing_im[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* pairing_im[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			fprintf(fp_angle, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n",w,s1kxky,s2kxky,x1kxky,x2kxky,p1kxky,p2kxky) ;
		}
		for(ie=1;ie<=Nener;ie++)
		{
			w = enstep*(ie-0.5) ;
			iee = Nener +ie + 1 ; 
			temp = 0;
			for(ix=0;ix<=nmom5;ix++)
			{
				if(temp<Akxky[iee][ix])
				{
					temp = Akxky[iee][ix];
					ix2 = ix ;
				}
			}
			kx = ix2*step2 ;
			ky = 1.0-(1.0-kx)*tan(M_PI*theta/180) ;
			ix1 = 0 ;
			iy1 = 0 ;
			count_x = 0 ;
			count_y = 0 ;

			for(ix2=0;ix2<=nmom-1;ix2++)
			{
				kx1 = ix2*step;
				kx3 = (ix2+1)*step;
				if((kx1<=fabs(kx))&&(fabs(kx)<=kx3)&&(count_x==0))
				{
					ix1 = ix2;
					count_x = count_x + 1 ;
				}
			}
			for(iy2=0;iy2<=nmom-1;iy2++)
			{
				ky1 = iy2*step;
				ky3 = (iy2+1)*step;
				if((ky1<=fabs(ky))&&(fabs(ky)<=ky3)&&(count_y==0))
				{
					iy1 = iy2;
					count_y = count_y + 1 ;
				}
			}
			ix3 = ix1+1;
			iy3 = iy1+1;
			kx1 = ix1*step;
			kx3 = (ix1+1)*step;
			ky1 = iy1*step;
			ky3 = (iy1+1)*step;

			s1kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* self_re[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* self_re[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* self_re[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* self_re[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			s2kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* self_im[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* self_im[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* self_im[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* self_im[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			x1kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* renor_re[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* renor_re[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* renor_re[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* renor_re[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			x2kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* renor_im[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* renor_im[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* renor_im[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* renor_im[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			p1kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* pairing_re[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* pairing_re[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* pairing_re[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* pairing_re[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			p2kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* pairing_im[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* pairing_im[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* pairing_im[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* pairing_im[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			fprintf(fp_angle, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n",w,s1kxky,s2kxky,x1kxky,x2kxky,p1kxky,p2kxky) ;
		}
		fclose(fp_angle) ;
	}



	sprintf(fid_angle2 , "sxp2");
	mkdir(fid_angle2,1017);
	for(itheta=0;itheta<=9;itheta++)
	{
		theta = itheta * 5;

		for(ix=0;ix<=nmom5;ix++)
		{
			kx = ix*step2 ;
			ky = 1.0-(1.0-kx)*tan(M_PI*theta/180) ;
			ix1 = 0 ;
			iy1 = 0 ;
			count_x = 0 ;
			count_y = 0 ;

			for(ix2=0;ix2<=nmom-1;ix2++)
			{
				kx1 = ix2*step;
				kx3 = (ix2+1)*step;
				if((kx1<=fabs(kx))&&(fabs(kx)<=kx3)&&(count_x==0))
				{
					ix1 = ix2;
					count_x = count_x + 1 ;
				}
			}

			if(count_x==0)
			{	
				printf("ix =%d\n",ix) ;
			}

			for(iy2=0;iy2<=nmom-1;iy2++)
			{
				ky1 = iy2*step;
				ky3 = (iy2+1)*step;
				if((ky1<=fabs(ky))&&(fabs(ky)<=ky3)&&(count_y==0))
				{
					iy1 = iy2;
					count_y = count_y + 1 ;
				}
			}
			if(count_y==0)
			{	
				printf("iy =%d\n",iy) ;
			}

			ix3 = ix1+1;
			iy3 = iy1+1;
			kx1 = ix1*step;
			kx3 = (ix1+1)*step;
			ky1 = iy1*step;
			ky3 = (iy1+1)*step;

			ie=1 ;
			w = -enstep*(ie-0.5) ;
			A11 = As[ie][ix1][iy1]-Ax[ie][ix1][iy1]  ;
			A31 = As[ie][ix3][iy1]-Ax[ie][ix3][iy1]  ;
			A13 = As[ie][ix1][iy3]-Ax[ie][ix1][iy3]  ;
			A33 = As[ie][ix3][iy3]-Ax[ie][ix3][iy3]  ;
			Akxky0[ix] = ((kx3-fabs(kx))*(ky3-fabs(ky))* A11+(fabs(kx)-kx1)*(ky3-fabs(ky))* A31+(kx3-fabs(kx))*(fabs(ky)-ky1)* A13+(fabs(kx)-kx1)*(fabs(ky)-ky1)* A33)/((kx3-kx1)*(kx3-kx1));
		}

		sprintf(fi_angle2, "sxp2/%d.dat",theta) ;
		fp_angle2 = fopen(fi_angle2,"w") ;

		temp = 0;
		for(ix=0;ix<=nmom5;ix++)
		{
			if(temp<Akxky0[ix])
			{
				temp = Akxky0[ix];
				ix2 = ix ;
			}
		}
		kx = ix2*step2 ;
		ky = 1.0-(1.0-kx)*tan(M_PI*theta/180) ;

		printf("theta = %d\tkx = %f\t ky = %f\n", theta,kx,ky);	
		fprintf(fp,"theta = %d\tkx = %f\t ky = %f\n", theta,kx,ky);	

		ix1 = 0 ;
		iy1 = 0 ;
		count_x = 0 ;
		count_y = 0 ;
		for(ix2=0;ix2<=nmom-1;ix2++)
		{
			kx1 = ix2*step;
			kx3 = (ix2+1)*step;
			if((kx1<=fabs(kx))&&(fabs(kx)<=kx3)&&(count_x==0))
			{
				ix1 = ix2;
				count_x = count_x + 1 ;
			}
		}
		for(iy2=0;iy2<=nmom-1;iy2++)
		{
			ky1 = iy2*step;
			ky3 = (iy2+1)*step;
			if((ky1<=fabs(ky))&&(fabs(ky)<=ky3)&&(count_y==0))
			{
				iy1 = iy2;
				count_y = count_y + 1 ;
			}
		}
		ix3 = ix1+1;
		iy3 = iy1+1;
		kx1 = ix1*step;
		kx3 = (ix1+1)*step;
		ky1 = iy1*step;
		ky3 = (iy1+1)*step;

		for(ie=1;ie<=Nener;ie++)
		{
			w = enstep*(ie-0.5) ;
			iee = Nener +ie + 1 ; 
			s1kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* self_re[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* self_re[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* self_re[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* self_re[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			s2kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* self_im[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* self_im[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* self_im[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* self_im[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			x1kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* renor_re[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* renor_re[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* renor_re[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* renor_re[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			x2kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* renor_im[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* renor_im[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* renor_im[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* renor_im[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			p1kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* pairing_re[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* pairing_re[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* pairing_re[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* pairing_re[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			p2kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* pairing_im[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* pairing_im[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* pairing_im[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* pairing_im[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			fprintf(fp_angle2, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n",w,s1kxky,s2kxky,x1kxky,x2kxky,p1kxky,p2kxky) ;
		}
		fclose(fp_angle2) ;

		sprintf(fi_angle3, "sxp2/imp%d.dat",theta) ;
		fp_angle3 = fopen(fi_angle3,"w") ;
		for(ie=1;ie<=Nener;ie++)
		{
			w = enstep*(ie-0.5) ;
			iee = Nener +ie + 1 ; 
			s1kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* self_imp_re[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* self_imp_re[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* self_imp_re[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* self_imp_re[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			s2kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* self_imp_im[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* self_imp_im[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* self_imp_im[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* self_imp_im[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			x1kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* renor_imp_re[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* renor_imp_re[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* renor_imp_re[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* renor_imp_re[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			x2kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* renor_imp_im[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* renor_imp_im[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* renor_imp_im[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* renor_imp_im[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			p1kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* pairing_imp_re[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* pairing_imp_re[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* pairing_imp_re[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* pairing_imp_re[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			p2kxky = ((kx3-fabs(kx))*(ky3-fabs(ky))* pairing_imp_im[ie][ix1][iy1]+(fabs(kx)-kx1)*(ky3-fabs(ky))* pairing_imp_im[ie][ix3][iy1]+(kx3-fabs(kx))*(fabs(ky)-ky1)* pairing_imp_im[ie][ix1][iy3]+(fabs(kx)-kx1)*(fabs(ky)-ky1)* pairing_imp_im[ie][ix3][iy3])/((kx3-kx1)*(kx3-kx1));
			fprintf(fp_angle3, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n",w,s1kxky,s2kxky,x1kxky,x2kxky,p1kxky,p2kxky) ;
		}
		fclose(fp_angle3) ;
	}


	
//	time
	time(&time_now);
	printf(ctime(&time_now));
	fprintf(fp,ctime(&time_now));
	running_time = time_now - time_0   ;
	printf("h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	
	fprintf(fp,"h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	


	fpp00 = fopen("A_nodal.dat","w");
	for(ix=0;ix<=nmom;ix++)
	{
//		if(ix%2==0)
//		{
			kx = step*(ix) ;
			iy = ix ;
			ky = step*(iy) ;
			for(ie=400;ie>=1;ie--)
			{
				w = -enstep*(ie-0.5) ;
				den0 = As[ie][ix][iy] - Ax[ie][ix][iy] ;
				fprintf(fpp00,"%f\t%f\t%f\n" ,kx,w,den0) ;
			}
			for(ie=1;ie<=400;ie++)
			{
				w = enstep*(ie-0.5) ;
				den0 = As[ie][ix][iy] + Ax[ie][ix][iy] ;
				fprintf(fpp00,"%f\t%f\t%f\n" ,kx,w,den0) ;
			}
//		}
	}

	fclose(fpp00) ;


	fpp = fopen("Aa.dat","w");
	iy = nmom ;
	ky = step*(iy) ;
	for(ix=0;ix<=nmom;ix++)
	{
//		if(ix%2==0)
//		{
			kx = step*(ix) ;
			for(ie=400;ie>=1;ie--)
			{
				w = -enstep*(ie-0.5) ;
				den0 = As[ie][ix][iy] - Ax[ie][ix][iy] ;
				fprintf(fpp,"%f\t%f\t%f\n" ,kx,w,den0) ;
			}
			for(ie=1;ie<=400;ie++)
			{
				w = enstep*(ie-0.5) ;
				den0 = As[ie][ix][iy] + Ax[ie][ix][iy] ;
				fprintf(fpp,"%f\t%f\t%f\n" ,kx,w,den0) ;
			}
//		}
	}

	fclose(fpp) ;


	fpp0 = fopen("A0.dat","w");
	iy = 0 ;
	ky = step*(iy) ;
	for(ix=0;ix<=nmom;ix++)
	{
//		if(ix%2==0)
//		{
			kx = step*(ix) ;
			for(ie=400;ie>=1;ie--)
			{
				w = -enstep*(ie-0.5) ;
				den0 = As[ie][ix][iy] - Ax[ie][ix][iy] ;
				fprintf(fpp0,"%f\t%f\t%f\n" ,kx,w,den0) ;
			}
			for(ie=1;ie<=400;ie++)
			{
				w = enstep*(ie-0.5) ;
				den0 = As[ie][ix][iy] + Ax[ie][ix][iy] ;
				fprintf(fpp0,"%f\t%f\t%f\n" ,kx,w,den0) ;
			}

//		}
	}

	fclose(fpp0) ;


// fourier den start
	for(ie=1;ie<=Nener;ie++) 
	{
		for(ix=1;ix<=2*nmom;ix++)
		{
			for(iy=1;iy<=2*nmom;iy++)	
			{
				if(ix<=nmom)
				{
					ix2 = ix - 1 ;
				}
				else
				{
					ix2 = 2*nmom - ix + 1 ;
				}
				if(iy<=nmom)
				{
					iy2 = iy - 1 ;
				}
				else
				{
					iy2 = 2*nmom - iy + 1 ;
				}
				As_temp[(2*ix-1)+(iy-1)*2*2*nmom] = As[ie][ix2][iy2] ;
				As_temp[(2*ix)+(iy-1)*2*2*nmom] = 0 ;
				Ax_temp[(2*ix-1)+(iy-1)*2*2*nmom] = Ax[ie][ix2][iy2] ;
				Ax_temp[(2*ix)+(iy-1)*2*2*nmom] = 0 ;
				Ap_temp[(2*ix-1)+(iy-1)*2*2*nmom] = Ap[ie][ix2][iy2] ;
				Ap_temp[(2*ix)+(iy-1)*2*2*nmom] = 0 ;
			}
		}
		fourn(As_temp,nmom2,2,-1) ;
		fourn(Ax_temp,nmom2,2,-1) ;
		fourn(Ap_temp,nmom2,2,-1) ;
		for(ix=0;ix<=nmom;ix++)
		{
			for(iy=0;iy<=nmom;iy++)	
			{
				ix2 = ix + 1 ;
				iy2 = iy + 1 ;
				As[ie][ix][iy] = As_temp[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
				Ax[ie][ix][iy] = Ax_temp[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
				Ap[ie][ix][iy] = Ap_temp[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
			}
		}
	}
// fourier den end

	for(ix=0;ix<=nmom;ix++)
	{
		kx = step*(ix) ;
		for(iy=0;iy<=ix;iy++)	
		{
			ky = step*(iy) ;
			#pragma omp parallel for private(ie)
			for(ie=1;ie<=Nener;ie++) 
			{
				w = enstep*(ie-0.5) ;
				s_sep_im[ie] = 0 ;
				x_sep_im[ie] = 0 ;
				p_sep_im[ie] = 0 ;

				for(je=1;je<=Nener;je++)
				{
					ener = enstep*(je-0.5) ;
					if((ie-je>=0)&&(ie+iwb<=Nener))
					{
						s_sep_im[ie] = s_sep_im[ie] - enstep*As[je][ix][iy]*((fermi_d(ener)-fermi_d(w+ener))*As[ie+je][ix][iy] 
																	+(fermi_d(-ener)-fermi_d(w-ener))*As[ie-je][ix][iy] ) ;
						x_sep_im[ie] = x_sep_im[ie] - enstep*Ax[je][ix][iy]*((fermi_d(ener)-fermi_d(w+ener))*Ax[ie+je][ix][iy] 
																	-(fermi_d(-ener)-fermi_d(w-ener))*Ax[ie-je][ix][iy] ) ;
						p_sep_im[ie] = p_sep_im[ie] - enstep*Ap[je][ix][iy]*((fermi_d(ener)-fermi_d(w+ener))*Ap[ie+je][ix][iy] 
																	-(fermi_d(-ener)-fermi_d(w-ener))*Ap[ie-je][ix][iy] ) ;
					}
					else if((ie-je<0)&&(ie+je<=Nener))
					{
						s_sep_im[ie] = s_sep_im[ie] - enstep*As[je][ix][iy]*((fermi_d(ener)-fermi_d(w+ener))*As[ie+je][ix][iy] 
																	+(fermi_d(-ener)-fermi_d(w-ener))*As[je-ie][ix][iy] ) ;
						x_sep_im[ie] = x_sep_im[ie] - enstep*Ax[je][ix][iy]*((fermi_d(ener)-fermi_d(w+ener))*Ax[ie+je][ix][iy] 
																	+(fermi_d(-ener)-fermi_d(w-ener))*Ax[je-ie][ix][iy] ) ;
						p_sep_im[ie] = p_sep_im[ie] - enstep*Ap[je][ix][iy]*((fermi_d(ener)-fermi_d(w+ener))*Ap[ie+je][ix][iy] 
																	+(fermi_d(-ener)-fermi_d(w-ener))*Ap[je-ie][ix][iy] ) ;
					}
					else if((ie-je>=0)&&(ie+je>Nener))
					{
						s_sep_im[ie] = s_sep_im[ie] - enstep*As[je][ix][iy]*((fermi_d(ener)-fermi_d(w+ener))*As[Nener][ix][iy] 
																	+(fermi_d(-ener)-fermi_d(w-ener))*As[ie-je][ix][iy] ) ;
						x_sep_im[ie] = x_sep_im[ie] - enstep*Ax[je][ix][iy]*((fermi_d(ener)-fermi_d(w+ener))*Ax[Nener][ix][iy] 
																	-(fermi_d(-ener)-fermi_d(w-ener))*Ax[ie-je][ix][iy] ) ;
						p_sep_im[ie] = p_sep_im[ie] - enstep*Ap[je][ix][iy]*((fermi_d(ener)-fermi_d(w+ener))*Ap[Nener][ix][iy] 
																	-(fermi_d(-ener)-fermi_d(w-ener))*Ap[ie-je][ix][iy] ) ;

					}
					else if((ie-je<0)&&(ie+je>Nener))
					{
						s_sep_im[ie] = s_sep_im[ie] - enstep*As[je][ix][iy]*((fermi_d(ener)-fermi_d(w+ener))*As[Nener][ix][iy] 
																	+(fermi_d(-ener)-fermi_d(w-ener))*As[je-ie][ix][iy] ) ;
						x_sep_im[ie] = x_sep_im[ie] - enstep*Ax[je][ix][iy]*((fermi_d(ener)-fermi_d(w+ener))*Ax[Nener][ix][iy] 
																	+(fermi_d(-ener)-fermi_d(w-ener))*Ax[je-ie][ix][iy] ) ;
						p_sep_im[ie] = p_sep_im[ie] - enstep*Ap[je][ix][iy]*((fermi_d(ener)-fermi_d(w+ener))*Ap[Nener][ix][iy] 
																	+(fermi_d(-ener)-fermi_d(w-ener))*Ap[je-ie][ix][iy] ) ;
					}
				}

				sep_s_im[ie][ix][iy] = M_PI*s_sep_im[ie] ;
				sep_s_im[ie][iy][ix] = M_PI*s_sep_im[ie] ;
				sep_x_im[ie][ix][iy] = M_PI*x_sep_im[ie] ;
				sep_x_im[ie][iy][ix] = M_PI*x_sep_im[ie] ;
				sep_p_im[ie][ix][iy] = M_PI*p_sep_im[ie] ;
				sep_p_im[ie][iy][ix] = M_PI*p_sep_im[ie] ;

			}

			#pragma omp parallel for private(ie)
			for(ie=1;ie<=Nener;ie++) 
			{
				s_sep_re[ie] = 0 ;
				x_sep_re[ie] = 0 ;
				p_sep_re[ie] = 0 ;
				w = enstep*(ie-0.5) ;
				for(je=1;je<=Nener;je++) 
				{
					ener = enstep*(je-0.5)  ;
					aa = (pow(ener,2)-pow(w,2))/(pow(pow(ener,2)-pow(w,2),2)+pow(deltkk,2)) ;
					s_sep_re[ie] = s_new_re[ie] + enstep*(1.0/M_PI) * (2*ener*sep_s_im[je][ix][iy]-2*w*sep_s_im[ie][ix][iy]) * aa;
					x_sep_re[ie] = x_new_re[ie] + enstep*(1.0/M_PI) * (2*ener*sep_x_im[je][ix][iy]-2*w*sep_x_im[ie][ix][iy]) * aa;
					p_sep_re[ie] = p_new_re[ie] + enstep*(1.0/M_PI) * (2*ener*sep_p_im[je][ix][iy]-2*w*sep_p_im[ie][ix][iy]) * aa;
				}

				sep_s_re[ie][ix][iy] = s_sep_re[ie] ;
				sep_s_re[ie][iy][ix] = s_sep_re[ie] ;
				sep_x_re[ie][ix][iy] = x_sep_re[ie];
				sep_x_re[ie][iy][ix] = x_sep_re[ie] ;
				sep_p_re[ie][ix][iy] = p_sep_re[ie] ;
				sep_p_re[ie][iy][ix] = p_sep_re[ie] ;
			}	

		}
	}	


// inverse fourier
	for(ie=1;ie<=Nener;ie++) 
	{
		for(ix=1;ix<=2*nmom;ix++)
		{
			for(iy=1;iy<=2*nmom;iy++)
			{
				if(ix<=nmom)
				{
					ix2 = ix - 1 ;
				}
				else
				{
					ix2 = 2*nmom - ix + 1 ;
				}
				if(iy<=nmom)
				{
					iy2 = iy - 1 ;
				}
				else
				{
					iy2 = 2*nmom - iy + 1 ;
				}
				sep_s_temp[(2*ix-1)+(iy-1)*2*2*nmom] = sep_s_re[ie][ix2][iy2] ;
				sep_s_temp[(2*ix)+(iy-1)*2*2*nmom] = sep_s_im[ie][ix2][iy2] ;
				sep_x_temp[(2*ix-1)+(iy-1)*2*2*nmom] = sep_x_re[ie][ix2][iy2] ;
				sep_x_temp[(2*ix)+(iy-1)*2*2*nmom] = sep_x_im[ie][ix2][iy2] ;
				sep_p_temp[(2*ix-1)+(iy-1)*2*2*nmom] = sep_p_re[ie][ix2][iy2] ;
				sep_p_temp[(2*ix)+(iy-1)*2*2*nmom] = sep_p_im[ie][ix2][iy2] ;
			}
		}
		fourn(sep_s_temp,nmom2,2,1) ;
		fourn(sep_x_temp,nmom2,2,1) ;
		fourn(sep_p_temp,nmom2,2,1) ;
		for(ix=0;ix<=nmom;ix++)
		{
			for(iy=0;iy<=nmom;iy++)	
			{
				ix2 = ix + 1 ;
				iy2 = iy + 1 ;
				sep_s_re[ie][ix][iy] = sep_s_temp[(2*ix2-1)+(iy2-1)*2*2*nmom]  ;
				sep_s_im[ie][ix][iy] = sep_s_temp[(2*ix2)+(iy2-1)*2*2*nmom] ;
				sep_x_re[ie][ix][iy] = sep_x_temp[(2*ix2-1)+(iy2-1)*2*2*nmom]  ;
				sep_x_im[ie][ix][iy] = sep_x_temp[(2*ix2)+(iy2-1)*2*2*nmom] ;
				sep_p_re[ie][ix][iy] = sep_p_temp[(2*ix2-1)+(iy2-1)*2*2*nmom]  ;
				sep_p_im[ie][ix][iy] = sep_p_temp[(2*ix2)+(iy2-1)*2*2*nmom] ;
			}
		}
	}






	fp_sep = fopen("sep_amp.dat","w");

	#pragma omp parallel for private(ie)
	for(ie=1;ie<=Nener;ie++)
	{
		w = (ie-0.5)*enstep ;
		s_sep_re[ie]=0 ; 
		s_sep_im[ie]=0 ; 
		x_sep_re[ie]=0 ; 
		x_sep_im[ie]=0 ; 
		p_sep_re[ie]=0 ; 
		p_sep_im[ie]=0 ; 
		for(ix=0;ix<=nmom;ix++)
		{
			kx = step*(ix) ;
			for(iy=0;iy<=nmom;iy++)	
			{
				ky = step*(iy) ;
				if((ix==0)&&(iy==0))
				{
					s_sep_re[ie] = s_sep_re[ie] + step*step*sep_s_re[ie][ix][iy] ;
					x_sep_re[ie] = x_sep_re[ie] + step*step*sep_x_re[ie][ix][iy] ;
					p_sep_re[ie] = p_sep_re[ie] + step*step*sep_p_re[ie][ix][iy] ;
					s_sep_im[ie] = s_sep_im[ie] + step*step*sep_s_im[ie][ix][iy] ;
					x_sep_im[ie] = x_sep_im[ie] + step*step*sep_x_im[ie][ix][iy] ;
					p_sep_im[ie] = p_sep_im[ie] + step*step*sep_p_im[ie][ix][iy] ;
				}
				else if((ix==nmom)&&(iy==nmom))
				{
					s_sep_re[ie] = s_sep_re[ie] + step*step*sep_s_re[ie][ix][iy] ;
					x_sep_re[ie] = x_sep_re[ie] + step*step*sep_x_re[ie][ix][iy] ;
					p_sep_re[ie] = p_sep_re[ie] + step*step*sep_p_re[ie][ix][iy] ;
					s_sep_im[ie] = s_sep_im[ie] + step*step*sep_s_im[ie][ix][iy] ;
					x_sep_im[ie] = x_sep_im[ie] + step*step*sep_x_im[ie][ix][iy] ;
					p_sep_im[ie] = p_sep_im[ie] + step*step*sep_p_im[ie][ix][iy] ;
				}
				else if((ix==0)&&(iy==nmom))
				{
					s_sep_re[ie] = s_sep_re[ie] + step*step*sep_s_re[ie][ix][iy] ;
					x_sep_re[ie] = x_sep_re[ie] + step*step*sep_x_re[ie][ix][iy] ;
					p_sep_re[ie] = p_sep_re[ie] + step*step*sep_p_re[ie][ix][iy] ;
					s_sep_im[ie] = s_sep_im[ie] + step*step*sep_s_im[ie][ix][iy] ;
					x_sep_im[ie] = x_sep_im[ie] + step*step*sep_x_im[ie][ix][iy] ;
					p_sep_im[ie] = p_sep_im[ie] + step*step*sep_p_im[ie][ix][iy] ;
				}
				else if((ix==nmom)&&(iy==0))
				{
					s_sep_re[ie] = s_sep_re[ie] + step*step*sep_s_re[ie][ix][iy] ;
					x_sep_re[ie] = x_sep_re[ie] + step*step*sep_x_re[ie][ix][iy] ;
					p_sep_re[ie] = p_sep_re[ie] + step*step*sep_p_re[ie][ix][iy] ;
					s_sep_im[ie] = s_sep_im[ie] + step*step*sep_s_im[ie][ix][iy] ;
					x_sep_im[ie] = x_sep_im[ie] + step*step*sep_x_im[ie][ix][iy] ;
					p_sep_im[ie] = p_sep_im[ie] + step*step*sep_p_im[ie][ix][iy] ;
				}
				else if((ix==0)&&(iy!=0))
				{
					s_sep_re[ie] = s_sep_re[ie] + 2*step*step*sep_s_re[ie][ix][iy] ;
					x_sep_re[ie] = x_sep_re[ie] + 2*step*step*sep_x_re[ie][ix][iy] ;
					p_sep_re[ie] = p_sep_re[ie] + 2*step*step*sep_p_re[ie][ix][iy] ;
					s_sep_im[ie] = s_sep_im[ie] + 2*step*step*sep_s_im[ie][ix][iy] ;
					x_sep_im[ie] = x_sep_im[ie] + 2*step*step*sep_x_im[ie][ix][iy] ;
					p_sep_im[ie] = p_sep_im[ie] + 2*step*step*sep_p_im[ie][ix][iy] ;
				}
				else if((ix!=0)&&(iy==0))
				{
					s_sep_re[ie] = s_sep_re[ie] + 2*step*step*sep_s_re[ie][ix][iy] ;
					x_sep_re[ie] = x_sep_re[ie] + 2*step*step*sep_x_re[ie][ix][iy] ;
					p_sep_re[ie] = p_sep_re[ie] + 2*step*step*sep_p_re[ie][ix][iy] ;
					s_sep_im[ie] = s_sep_im[ie] + 2*step*step*sep_s_im[ie][ix][iy] ;
					x_sep_im[ie] = x_sep_im[ie] + 2*step*step*sep_x_im[ie][ix][iy] ;
					p_sep_im[ie] = p_sep_im[ie] + 2*step*step*sep_p_im[ie][ix][iy] ;
				}
				else if((ix==nmom)&&(iy!=nmom))
				{
					s_sep_re[ie] = s_sep_re[ie] + 2*step*step*sep_s_re[ie][ix][iy] ;
					x_sep_re[ie] = x_sep_re[ie] + 2*step*step*sep_x_re[ie][ix][iy] ;
					p_sep_re[ie] = p_sep_re[ie] + 2*step*step*sep_p_re[ie][ix][iy] ;
					s_sep_im[ie] = s_sep_im[ie] + 2*step*step*sep_s_im[ie][ix][iy] ;
					x_sep_im[ie] = x_sep_im[ie] + 2*step*step*sep_x_im[ie][ix][iy] ;
					p_sep_im[ie] = p_sep_im[ie] + 2*step*step*sep_p_im[ie][ix][iy] ;
				}
				else if((ix!=nmom)&&(iy==nmom))
				{
					s_sep_re[ie] = s_sep_re[ie] + 2*step*step*sep_s_re[ie][ix][iy] ;
					x_sep_re[ie] = x_sep_re[ie] + 2*step*step*sep_x_re[ie][ix][iy] ;
					p_sep_re[ie] = p_sep_re[ie] + 2*step*step*sep_p_re[ie][ix][iy] ;
					s_sep_im[ie] = s_sep_im[ie] + 2*step*step*sep_s_im[ie][ix][iy] ;
					x_sep_im[ie] = x_sep_im[ie] + 2*step*step*sep_x_im[ie][ix][iy] ;
					p_sep_im[ie] = p_sep_im[ie] + 2*step*step*sep_p_im[ie][ix][iy] ;
				}
				else
				{
					s_sep_re[ie] = s_sep_re[ie] + 4*step*step*sep_s_re[ie][ix][iy] ;
					x_sep_re[ie] = x_sep_re[ie] + 4*step*step*sep_x_re[ie][ix][iy] ;
					p_sep_re[ie] = p_sep_re[ie] + 4*step*step*sep_p_re[ie][ix][iy] ;
					s_sep_im[ie] = s_sep_im[ie] + 4*step*step*sep_s_im[ie][ix][iy] ;
					x_sep_im[ie] = x_sep_im[ie] + 4*step*step*sep_x_im[ie][ix][iy] ;
					p_sep_im[ie] = p_sep_im[ie] + 4*step*step*sep_p_im[ie][ix][iy] ;
				}
			}
		}
	}
	
	for(ie=1;ie<=Nener;ie++)
	{
		w = (ie-0.5)*enstep ;
		fprintf(fp_sep,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",w,s_sep_re[ie],s_sep_im[ie],x_sep_re[ie],x_sep_im[ie],p_sep_re[ie],p_sep_im[ie]) ;
	}
	fclose(fp_sep);


	fp_sep_s = fopen("sep_s.dat","w") ;
	fp_sep_x = fopen("sep_x.dat","w") ;
	fp_sep_p = fopen("sep_p.dat","w") ;
	for(ie=1;ie<=500;ie++)
	{
		w = (ie-0.5)*enstep ;
		for(ix=0;ix<=nmom;ix++)
		{
			kx = step*(ix) ;
			iy = ix ;
			ky = step*(iy) ;
			fprintf(fp_sep_s, "%f\t%f\t%f\n",w,kx,sep_s_im[ie][ix][iy]) ;
			fprintf(fp_sep_x, "%f\t%f\t%f\n",w,kx,sep_x_im[ie][ix][iy]) ;
			fprintf(fp_sep_p, "%f\t%f\t%f\n",w,kx,sep_p_im[ie][ix][iy]) ;
		}
	}

	fclose(fp_sep_s);
	fclose(fp_sep_x);
	fclose(fp_sep_p);



	fclose(fp);	




}





double xi_band(double kx, double ky)
{
	double xi_b ;
	xi_b = -2*t1 *(cos(M_PI*kx) + cos(M_PI*ky)) + 4*t2 *cos(M_PI*kx) * cos(M_PI*ky)- 2*t3*(cos(2*M_PI*kx) + cos(2*M_PI*ky)) - chem ;
	return xi_b ;
}

double fermi_d(double ener)
{
	double fer ;

	if(fabs(ener)<0.0000001)
	{
		fer = 0.5 ;
	}
	else
	{
		fer = 1.0/(exp((ener)/Temp) + 1) ;
	}
	
	return fer ;
}

double bose_d(double wboson)
{
	double bos ;

	if(fabs(wboson)<0.0000001)
	{
		bos = 0.0 ;
	}
	else
	{
		bos = 1.0/(exp((wboson)/Temp) - 1) ;
	}
	
	return bos ;
}


double d_wave(double kx, double ky)
{
	return (cos(M_PI*kx)-cos(M_PI*ky))/2.0 ;
}



double Vimp(double qx, double qy)
{
	double Rq;

	Rq = 2.0*M_PI*kapa_imp*V0/pow((pow(M_PI*qx,2)+pow(M_PI*qy,2))+pow(kapa_imp,2),1.5) ;

	return Rq ;
}



























