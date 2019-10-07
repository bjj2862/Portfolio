#ifndef _cal_function_H
#define _cal_function_H

#include "some_function.h"

void Cal_Green_function(complex_3D G, double_3D A, complex_3D self_eff, complex_3D self_imp)
{
	int ie, ix,iy  ;

	double ener, kx, ky ;
	double *ws1, *xx1, *de_r, *de_i ;

	ws1 = dvector(0,Nener);
	xx1 =  dvector(0,Nener);
	de_r = dvector(0,Nener);
	de_i = dvector(0,Nener);


	for(ix=0;ix<=nmom;ix++)
	{
		kx = momentum(ix) ;
		for(iy=0;iy<=ix;iy++)	
		{
			ky = momentum(iy) ;
			#pragma omp parallel for private(ie)
			for(ie=1;ie<=Nener;ie++) 
			{

//				(n1+i*n2)/(de1-i*de2)

				ws1[ie] = energy_fermion(ie) - (self_eff.s1[ie][ix][iy] + self_imp.s1[ie][ix][iy]);
				xx1[ie] = xi_band(kx,ky) + (self_eff.x1[ie][ix][iy] + self_imp.x1[ie][ix][iy]) ;

				de_r[ie] = pow(ws1[ie],2) - pow(self_eff.s2[ie][ix][iy]+self_imp.s2[ie][ix][iy],2) - pow(xx1[ie],2) + pow(self_eff.x2[ie][ix][iy]+self_imp.x2[ie][ix][iy],2) - pow(self_eff.p1[ie][ix][iy]+self_imp.p1[ie][ix][iy],2) + pow(self_eff.p2[ie][ix][iy]+self_imp.p2[ie][ix][iy],2) ;
				de_i[ie] = 2*ws1[ie]*(self_eff.s2[ie][ix][iy]+self_imp.s2[ie][ix][iy]) + 2*xx1[ie]*(self_eff.x2[ie][ix][iy]+self_imp.x2[ie][ix][iy]) + 2*(self_eff.p1[ie][ix][iy]+self_imp.p1[ie][ix][iy])*(self_eff.p2[ie][ix][iy]+self_imp.p2[ie][ix][iy]) ;

				G.s1[ie][ix][iy] = (ws1[ie] * de_r[ie] + (self_eff.s2[ie][ix][iy]+self_imp.s2[ie][ix][iy])*de_i[ie])/(pow(de_r[ie],2)+pow(de_i[ie],2));
				G.x1[ie][ix][iy] = (xx1[ie] * de_r[ie] - (self_eff.x2[ie][ix][iy]+self_imp.x2[ie][ix][iy])*de_i[ie] )/(pow(de_r[ie],2)+pow(de_i[ie],2));
				G.p1[ie][ix][iy] = ((self_eff.p1[ie][ix][iy]+self_imp.p1[ie][ix][iy]) * de_r[ie] - (self_eff.p2[ie][ix][iy]+self_imp.p2[ie][ix][iy])*de_i[ie] )/(pow(de_r[ie],2)+pow(de_i[ie],2));

				G.s2[ie][ix][iy] = (ws1[ie] * de_i[ie] - (self_eff.s2[ie][ix][iy]+self_imp.s2[ie][ix][iy])*de_r[ie])/(pow(de_r[ie],2)+pow(de_i[ie],2));
				G.x2[ie][ix][iy] = (xx1[ie] * de_i[ie] + (self_eff.x2[ie][ix][iy]+self_imp.x2[ie][ix][iy])*de_r[ie] )/(pow(de_r[ie],2)+pow(de_i[ie],2));
				G.p2[ie][ix][iy] = ((self_eff.p1[ie][ix][iy]+self_imp.p1[ie][ix][iy]) * de_i[ie] + (self_eff.p2[ie][ix][iy]+self_imp.p2[ie][ix][iy])*de_r[ie] )/(pow(de_r[ie],2)+pow(de_i[ie],2));

				A.s[ie][ix][iy] = - (1.0/M_PI) * G.s2[ie][ix][iy];
				A.x[ie][ix][iy] = - (1.0/M_PI) * G.x2[ie][ix][iy];
				A.p[ie][ix][iy] = - (1.0/M_PI) * G.p2[ie][ix][iy];

			}
		}
	}

	free_dvector(ws1,0,Nener);
	free_dvector(xx1,0,Nener);
	free_dvector(de_r,0,Nener);
	free_dvector(de_i,0,Nener);
}

void cal_eff_self_energy(complex_3D self_eff, double_3D A, double ***Vq_eff)
{
	int  ix,iy,ie,iwb, je ;
	double kx, ky, w, wboson, ener, aa;
	double *s_new, *x_new, *p_new;
//	int Nthread=omp_get_thread_num();
//	int ener_f(Nthread), ener_b(Nthread);


	s_new = dvector(0,Nener);
	x_new = dvector(0,Nener);
	p_new = dvector(0,Nener);

	for(ix=0;ix<=nmom;ix++)
	{
		kx = momentum(ix) ;
		for(iy=0;iy<=ix;iy++)	
		{
			ky = momentum(iy) ;

			#pragma omp parallel for private(ie,iwb,w,wboson)
			for(ie=1;ie<=Nener;ie++) 
			{
				w = energy_fermion(ie) ;
				s_new[ie] = 0 ;
				x_new[ie] = 0 ;
				p_new[ie] = 0 ;

				for(iwb=iwbs;iwb<=iwbc;iwb++)
				{
					wboson = energy_boson(iwb) ;

					if((ie-iwb>=0)&&(ie+iwb<=Nener))
					{
						s_new[ie] = s_new[ie]+enstep*Vq_eff[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*A.s[ie-iwb][ix][iy] 
											-(fermi_d(w+wboson)+bose_d(+wboson))*A.s[ie+iwb][ix][iy]) ;	
						x_new[ie] = x_new[ie]+enstep*Vq_eff[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*A.x[ie-iwb][ix][iy] 
											-(fermi_d(w+wboson)+bose_d(+wboson))*A.x[ie+iwb][ix][iy]) ;	
						p_new[ie] = p_new[ie]+enstep*Vq_eff[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*A.p[ie-iwb][ix][iy] 
											-(fermi_d(w+wboson)+bose_d(+wboson))*A.p[ie+iwb][ix][iy]) ;	
  					}
	  				else if((ie-iwb<0)&&(ie+iwb<=Nener))
  					{
						s_new[ie] = s_new[ie]+enstep*Vq_eff[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*A.s[iwb-ie][ix][iy] 
											-(fermi_d(w+wboson)+bose_d(+wboson))*A.s[ie+iwb][ix][iy]) ;	
						x_new[ie] = x_new[ie]+enstep*Vq_eff[iwb][ix][iy]*(-(fermi_d(w-wboson)+bose_d(-wboson))*A.x[iwb-ie][ix][iy] 
											-(fermi_d(w+wboson)+bose_d(+wboson))*A.x[ie+iwb][ix][iy]) ;	
						p_new[ie] = p_new[ie]+enstep*Vq_eff[iwb][ix][iy]*(-(fermi_d(w-wboson)+bose_d(-wboson))*A.p[iwb-ie][ix][iy] 
										-(fermi_d(w+wboson)+bose_d(+wboson))*A.p[ie+iwb][ix][iy]) ;	
  					}
  					else if((ie-iwb>=0)&&(ie+iwb>Nener))
  					{
						s_new[ie] = s_new[ie]+enstep*Vq_eff[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*A.s[ie-iwb][ix][iy] 
											-(fermi_d(w+wboson)+bose_d(+wboson))*A.s[Nener][ix][iy]) ;
						x_new[ie] = x_new[ie]+enstep*Vq_eff[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*A.x[ie-iwb][ix][iy] 
											-(fermi_d(w+wboson)+bose_d(+wboson))*A.x[Nener][ix][iy]) ;
						p_new[ie] = p_new[ie]+enstep*Vq_eff[iwb][ix][iy]*((fermi_d(w-wboson)+bose_d(-wboson))*A.p[ie-iwb][ix][iy] 
											-(fermi_d(w+wboson)+bose_d(+wboson))*A.p[Nener][ix][iy]) ;	
					}

				}

				self_eff.s2[ie][ix][iy] = M_PI*s_new[ie] ;
				self_eff.x2[ie][ix][iy] = M_PI*x_new[ie];
				self_eff.p2[ie][ix][iy] = M_PI*p_new[ie] ;
			}

			#pragma omp parallel for private(ie,je,w,ener,aa)
			for(ie=1;ie<=Nener;ie++) 
			{
				s_new[ie] = 0 ;
				x_new[ie] = 0 ;
				p_new[ie] = 0 ;
				w = energy_fermion(ie) ;
				for(je=1;je<=Nener;je++) 
				{
					ener = energy_fermion(je)  ;
					aa = (pow(ener,2)-pow(w,2))/(pow(pow(ener,2)-pow(w,2),2)+pow(deltkk,2)) ;
					s_new[ie] = s_new[ie] + enstep*(1.0/M_PI) * (2*w*self_eff.s2[je][ix][iy]-2*w*self_eff.s2[ie][ix][iy]) * aa;
					x_new[ie] = x_new[ie] + enstep*(1.0/M_PI) * (2*ener*self_eff.x2[je][ix][iy]-2*w*self_eff.x2[ie][ix][iy]) * aa;
					p_new[ie] = p_new[ie] + enstep*(1.0/M_PI) * (2*ener*self_eff.p2[je][ix][iy]-2*w*self_eff.p2[ie][ix][iy]) * aa;
				}
				self_eff.s1[ie][ix][iy] = s_new[ie] ;
				self_eff.x1[ie][ix][iy] = x_new[ie];
				self_eff.p1[ie][ix][iy] = p_new[ie] ;
			}

		}
	}	

	free_dvector(s_new,0,Nener);
	free_dvector(x_new,0,Nener);
	free_dvector(p_new,0,Nener);

}


void cal_imp_self_energy(complex_3D self_imp, complex_3D G, double **Vq_imp)
{
	int  ix,iy,ie ;
	double kx, ky, w;

	for(ix=0;ix<=nmom;ix++)
	{
		kx = momentum(ix) ;
		for(iy=0;iy<=ix;iy++)	
		{
			ky = momentum(iy) ;

			#pragma omp parallel for private(ie,w)
			for(ie=1;ie<=Nener;ie++) 
			{
				w = energy_fermion(ie) ;

				self_imp.s1[ie][ix][iy] = Vq_imp[ix][iy]*G.s1[ie][ix][iy] ;
				self_imp.x1[ie][ix][iy] = Vq_imp[ix][iy]*G.x1[ie][ix][iy] ;
				self_imp.p1[ie][ix][iy] = -Vq_imp[ix][iy]*G.p1[ie][ix][iy] ;

				self_imp.s2[ie][ix][iy] = Vq_imp[ix][iy]*G.s2[ie][ix][iy] ;
				self_imp.x2[ie][ix][iy] = Vq_imp[ix][iy]*G.x2[ie][ix][iy] ;
				self_imp.p2[ie][ix][iy] = -Vq_imp[ix][iy]*G.p2[ie][ix][iy] ;
			}
		}
	}	


}


complex_1D Cal_DOS(complex_3D G, double_3D A)
{
	int ie, ix,iy  ;

	double ener, kx, ky ;
	complex_1D DOS;
	DOS.s1 = dvector(0,Nener) ;
	DOS.s2 = dvector(0,Nener) ;
	DOS.x1 = dvector(0,Nener) ;
	DOS.x2 = dvector(0,Nener) ;

//       ppppp
	#pragma omp parallel for private(ie,ener,ix,iy,kx,ky)
	for(ie=1;ie<=Nener;ie++) 
	{
		ener = energy_fermion(ie) ;
		DOS.s1[ie] = 0 ;
		DOS.x1[ie] = 0 ;
		DOS.s2[ie] = 0 ;
		DOS.x2[ie] = 0 ;
	
		for(ix=0;ix<=nmom;ix++)
		{
			kx = momentum(ix) ;
			for(iy=0;iy<=ix;iy++)	
			{
				ky = momentum(iy) ;
				if((ix==0)&&(iy==0))
				{
					DOS.s1[ie] = DOS.s1[ie] + 1*step*step*A.s[ie][ix][iy]/2.0 ;
					DOS.x1[ie] = DOS.x1[ie] + 1*step*step*A.x[ie][ix][iy]/2.0 ;
					DOS.s2[ie] = DOS.s2[ie] - 1*M_PI*step*step*G.s1[ie][ix][iy]/2.0 ;
					DOS.x2[ie] = DOS.x2[ie] - 1*M_PI*step*step*G.x1[ie][ix][iy]/2.0 ;
				}
				else if((ix==nmom)&&(iy==nmom))
				{
					DOS.s1[ie] = DOS.s1[ie] + 1*step*step*A.s[ie][ix][iy]/2.0 ;
					DOS.x1[ie] = DOS.x1[ie] + 1*step*step*A.x[ie][ix][iy]/2.0 ;
					DOS.s2[ie] = DOS.s2[ie] - 1*M_PI*step*step*G.s1[ie][ix][iy]/2.0 ;
					DOS.x2[ie] = DOS.x2[ie] - 1*M_PI*step*step*G.x1[ie][ix][iy]/2.0 ;
				}
				else if((ix==nmom)&&(iy==0))
				{
					DOS.s1[ie] = DOS.s1[ie] + 2*step*step*A.s[ie][ix][iy]/2.0 ;
					DOS.x1[ie] = DOS.x1[ie] + 2*step*step*A.x[ie][ix][iy]/2.0 ;
					DOS.s2[ie] = DOS.s2[ie] - 2*M_PI*step*step*G.s1[ie][ix][iy]/2.0 ;
					DOS.x2[ie] = DOS.x2[ie] - 2*M_PI*step*step*G.x1[ie][ix][iy]/2.0 ;
				}
				else if((ix!=0)&&(iy==0))
				{
					DOS.s1[ie] = DOS.s1[ie] + 4*step*step*A.s[ie][ix][iy]/2.0 ;
					DOS.x1[ie] = DOS.x1[ie] + 4*step*step*A.x[ie][ix][iy]/2.0 ;
					DOS.s2[ie] = DOS.s2[ie] - 4*M_PI*step*step*G.s1[ie][ix][iy]/2.0 ;
					DOS.x2[ie] = DOS.x2[ie] - 4*M_PI*step*step*G.x1[ie][ix][iy]/2.0 ;
				}
				else if((ix==nmom)&&(iy!=nmom))
				{
					DOS.s1[ie] = DOS.s1[ie] + 4*step*step*A.s[ie][ix][iy]/2.0 ;
					DOS.x1[ie] = DOS.x1[ie] + 4*step*step*A.x[ie][ix][iy]/2.0 ;
					DOS.s2[ie] = DOS.s2[ie] - 4*M_PI*step*step*G.s1[ie][ix][iy]/2.0 ;
					DOS.x2[ie] = DOS.x2[ie] - 4*M_PI*step*step*G.x1[ie][ix][iy]/2.0 ;
				}
				else if(ix==iy)
				{
					DOS.s1[ie] = DOS.s1[ie] + 4*step*step*A.s[ie][ix][iy]/2.0 ;
					DOS.x1[ie] = DOS.x1[ie] + 4*step*step*A.x[ie][ix][iy]/2.0 ;
					DOS.s2[ie] = DOS.s2[ie] - 4*M_PI*step*step*G.s1[ie][ix][iy]/2.0 ;
					DOS.x2[ie] = DOS.x2[ie] - 4*M_PI*step*step*G.x1[ie][ix][iy]/2.0 ;
				}

				else
				{
					DOS.s1[ie] = DOS.s1[ie] + 8*step*step*A.s[ie][ix][iy]/2.0 ;
					DOS.x1[ie] = DOS.x1[ie] + 8*step*step*A.x[ie][ix][iy]/2.0 ;
					DOS.s2[ie] = DOS.s2[ie] - 8*M_PI*step*step*G.s1[ie][ix][iy]/2.0 ;
					DOS.x2[ie] = DOS.x2[ie] - 8*M_PI*step*step*G.x1[ie][ix][iy]/2.0 ;
				}
			}
		}
	}

	return DOS ;

}


#endif 

