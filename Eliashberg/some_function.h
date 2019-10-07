#ifndef _some_function_H
#define _some_function_H

double momentum(int i)
{
	return step*i;
}
double energy_fermion(int i)
{
	return enstep*(i-0.5);
}
double energy_boson(int i)
{
	return enstep*i;
}

int four_mom(int i)
{
	int j ;
	if(i<=nmom)
	{
		j = i - 1 ;
	}
	else
	{
		j = 2*nmom - i + 1 ;
	}
	return j ;
}
void fourn_m_to_r(double **AA, int symm)
{
	FILE *fp ;
	int ix, iy, ix2, iy2;
	int i;

	unsigned long *nmom2 ;
	float *temp_vec;
	nmom2 = lvector(1,2);
	nmom2[1] = 2*nmom ;
    nmom2[2] = 2*nmom ;	

	temp_vec = vector(1,2*4*nmom*nmom);

	for(ix=1;ix<=2*nmom;ix++)
	{
		ix2 = four_mom(ix);
		for(iy=1;iy<=2*nmom;iy++)	
		{
			iy2 = four_mom(iy);
			if(iy2<=ix2)
			{
				temp_vec[(2*ix-1)+(iy-1)*2*2*nmom]= AA[ix2][iy2] ;
			}
			else 
			{
				if(symm==1)
				{
						temp_vec[(2*ix-1)+(iy-1)*2*2*nmom]= AA[iy2][ix2] ; 
				}
				else if(symm==-1)
				{
						temp_vec[(2*ix-1)+(iy-1)*2*2*nmom]= -AA[iy2][ix2] ; 
				}
			}
			temp_vec[2*ix+(iy-1)*2*2*nmom] = 0 ;

		}
	}
	fourn(temp_vec,nmom2,2,-1) ;
	for(ix=0;ix<=nmom;ix++)
	{
		for(iy=0;iy<=ix;iy++)	
		{
			ix2 = ix + 1 ;
			iy2 = iy + 1 ;
			AA[ix][iy] = temp_vec[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
		}
	}

	free_lvector(nmom2,1,2);
	free_vector(temp_vec,1,2*4*nmom*nmom) ;

}


void fourn_m_to_r2(double **G1, double **G2,int symm)
{
	int ix, iy, ix2, iy2;
	
	unsigned long *nmom2 ;
	float *temp_vec;
	nmom2 = lvector(1,2);
	nmom2[1] = 2*nmom ;
    nmom2[2] = 2*nmom ;	

	temp_vec = vector(1,2*4*nmom*nmom);

	for(ix=1;ix<=2*nmom;ix++)
	{
		ix2 = four_mom(ix);
		for(iy=1;iy<=2*nmom;iy++)	
		{
			iy2 = four_mom(iy);
			if(iy2<=ix2)
			{
				temp_vec[(2*ix-1)+(iy-1)*2*2*nmom]= G1[ix2][iy2] ; 
				temp_vec[2*ix+(iy-1)*2*2*nmom] = G2[ix2][iy2] ;
			}
			else
			{
				if(symm==1)
				{
					temp_vec[(2*ix-1)+(iy-1)*2*2*nmom]= G1[iy2][ix2] ; 
					temp_vec[2*ix+(iy-1)*2*2*nmom] = G2[iy2][ix2] ;
				}
				else if(symm==-1)
				{
					temp_vec[(2*ix-1)+(iy-1)*2*2*nmom]= -G1[iy2][ix2] ; 
					temp_vec[2*ix+(iy-1)*2*2*nmom] = -G2[iy2][ix2] ;
				}

			}
		}
	}
	fourn(temp_vec,nmom2,2,-1) ;
	for(ix=0;ix<=nmom;ix++)
	{
		for(iy=0;iy<=ix;iy++)	
		{
			ix2 = ix + 1 ;
			iy2 = iy + 1 ;
			G1[ix][iy] = temp_vec[(2*ix2-1)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
			G2[ix][iy] = temp_vec[(2*ix2)+(iy2-1)*2*2*nmom] / pow(2*nmom,2) ;
		}
	}

	free_lvector(nmom2,1,2);
	free_vector(temp_vec,1,2*4*nmom*nmom) ;

}


void fourn_r_to_m(double **s1, double **s2, int symm)
{
	int ix, iy, ix2, iy2;
	
	unsigned long *nmom2 ;
	float *temp_vec;
	nmom2 = lvector(1,2);
	nmom2[1] = 2*nmom ;
    nmom2[2] = 2*nmom ;	

	temp_vec = vector(1,2*4*nmom*nmom);

	for(ix=1;ix<=2*nmom;ix++)
	{
		ix2 = four_mom(ix);
		for(iy=1;iy<=2*nmom;iy++)	
		{
			iy2 = four_mom(iy);

			if(iy2<=ix2)
			{
				temp_vec[(2*ix-1)+(iy-1)*2*2*nmom]= s1[ix2][iy2] ; 
				temp_vec[2*ix+(iy-1)*2*2*nmom] = s2[ix2][iy2] ;
			}
			else
			{
				if(symm==1)
				{
					temp_vec[(2*ix-1)+(iy-1)*2*2*nmom]= s1[iy2][ix2] ; 
					temp_vec[2*ix+(iy-1)*2*2*nmom] = s2[iy2][ix2] ;
				}
				else if(symm==-1)
				{
					temp_vec[(2*ix-1)+(iy-1)*2*2*nmom]= -s1[iy2][ix2] ; 
					temp_vec[2*ix+(iy-1)*2*2*nmom] = -s2[iy2][ix2] ;
				}

			}
		}
	}

	fourn(temp_vec,nmom2,2,1) ;
	for(ix=0;ix<=nmom;ix++)
	{
		for(iy=0;iy<=ix;iy++)	
		{
			ix2 = ix + 1 ;
			iy2 = iy + 1 ;
			s1[ix][iy] = temp_vec[(2*ix2-1)+(iy2-1)*2*2*nmom] ;
			s2[ix][iy] = temp_vec[(2*ix2)+(iy2-1)*2*2*nmom] ;

		}
	}

	free_lvector(nmom2,1,2);
	free_vector(temp_vec,1,2*4*nmom*nmom) ;

}
double *cal_difference(complex_3D self, complex_3D self_new)
{
	int i, ix, iy, ie ;
	double *difference ;
	double temp_s, temp_x, temp_p;

	for(ix=0;ix<=nmom;ix++)
	{
		for(iy=0;iy<=ix;iy++)	
		{
			for(ie=1;ie<=Nener;ie++) 
			{
				temp_s += pow(self_new.s1[ie][ix][iy]-self.s1[ie][ix][iy],2) + pow(self_new.s2[ie][ix][iy]-self.s2[ie][ix][iy],2) ;
				temp_x += pow(self_new.x1[ie][ix][iy]-self.x1[ie][ix][iy],2) + pow(self_new.x2[ie][ix][iy]-self.x2[ie][ix][iy],2) ;
				temp_p += pow(self_new.p1[ie][ix][iy]-self.p1[ie][ix][iy],2) + pow(self_new.p2[ie][ix][iy]-self.p2[ie][ix][iy],2) ;
			}
		}
	}
		
	difference = dvector(0,3);
	difference[1] = 2*sqrt(temp_s) * pow(step,2) * enstep ;
	difference[2] = 2*sqrt(temp_x) * pow(step,2) * enstep ;
	difference[3] = 2*sqrt(temp_p) * pow(step,2) * enstep ;
	difference[0] = difference[1]+difference[2]+difference[3];
	return difference;
}
void mixing(complex_3D self, complex_3D self_new, double mixing_rate)
{
	int ix, iy, ie ;
	for(ix=0;ix<=nmom;ix++)
	{
		for(iy=0;iy<=ix;iy++)	
		{
			for(ie=1;ie<=Nener;ie++) 
			{
					self.s1[ie][ix][iy] = (1.0-mixing_rate)*self.s1[ie][ix][iy] + mixing_rate*self_new.s1[ie][ix][iy] ;
					self.s2[ie][ix][iy] = (1.0-mixing_rate)*self.s2[ie][ix][iy] + mixing_rate*self_new.s2[ie][ix][iy] ;
					self.x1[ie][ix][iy] = (1.0-mixing_rate)*self.x1[ie][ix][iy] + mixing_rate*self_new.x1[ie][ix][iy] ;
					self.x2[ie][ix][iy] = (1.0-mixing_rate)*self.x2[ie][ix][iy] + mixing_rate*self_new.x2[ie][ix][iy] ;
					self.p1[ie][ix][iy] = (1.0-mixing_rate)*self.p1[ie][ix][iy] + mixing_rate*self_new.p1[ie][ix][iy] ;
					self.p2[ie][ix][iy] = (1.0-mixing_rate)*self.p2[ie][ix][iy] + mixing_rate*self_new.p2[ie][ix][iy] ;
			}	
		}
	}		
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

#endif


