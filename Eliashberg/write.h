
#ifndef _write_H
#define _write_H

void write_DOS(complex_1D DOS)
{
	FILE *fp_DOS;
	int ie;
	double ener;

	fp_DOS = fopen("DOS.d","w");

	for(ie=Nener;ie>=1;ie--) 
	{	
		ener = -energy_fermion(ie) ;
		fprintf(fp_DOS,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",ener,DOS.s1[ie],-DOS.s2[ie],-DOS.x1[ie],DOS.x2[ie],(DOS.s1[ie] - DOS.x1[ie]),(-DOS.s2[ie] + DOS.x2[ie])) ;
	}
	for(ie=1;ie<=Nener;ie++) 
	{	
		ener = energy_fermion(ie) ;
		fprintf(fp_DOS,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",ener,DOS.s1[ie],DOS.s2[ie],DOS.x1[ie],DOS.x2[ie],(DOS.s1[ie] + DOS.x1[ie]),(+DOS.s2[ie] + DOS.x2[ie])) ;
	}

	fclose(fp_DOS) ;
}

void write_minus_spectral_function(FILE *fp, double_3D A)
{
	FILE *fpmi ;
	int ie, ix, iy, Ami ;
	double kx, ky, ener;

	fpmi = fopen("Aminus.d","w") ;
	Ami = 0 ;
	for(ie=1;ie<=Nener;ie++) 
	{
		ener = energy_fermion(ie) ;
		for(ix=0;ix<=nmom;ix++)
		{
			kx = momentum(ix) ;
			for(iy=0;iy<=ix;iy++)	
			{
				ky = momentum(iy) ;
				if(A.s[ie][ix][iy]+A.x[ie][ix][iy]<0)
				{
					fprintf(fpmi,"%f\t%f\t%f\t%f\n",kx,ky,ener,A.s[ie][ix][iy]+A.x[ie][ix][iy]) ;
					Ami = Ami + 1 ;
				}
				if(A.s[ie][ix][iy]-A.x[ie][ix][iy]<0)
				{
					fprintf(fpmi,"%f\t%f\t%f\t%f\n",kx,ky,-ener,A.s[ie][ix][iy]-A.x[ie][ix][iy]) ;
					Ami = Ami + 1 ;
				}
			}
		}
	}
	printf("number of the minus A : %d\n", Ami) ;
	fprintf(fp,"number of the minus A : %d\n", Ami) ;
	fclose(fpmi) ;

}


void write_Fermi_surface(double_3D A)
{
	FILE *fp_fs ;
	int ix, iy, ie=1 ;
	double kx, ky, den0 ;

	fp_fs = fopen("FS.dat","w");
	for(ix=0;ix<=nmom;ix++)	
	{
		kx = momentum(ix) ;
			
		for(iy=0;iy<=nmom;iy++)	
		{
			ky = momentum(iy) ;
			if(iy<=ix)
			{
				den0 = A.s[ie][ix][iy] - A.x[ie][ix][iy] ;
			}
			else
			{
				den0 = A.s[ie][iy][ix] - A.x[ie][iy][ix] ;
			}
			fprintf(fp_fs,"%f\t%f\t%f\n", kx,ky,den0) ;
		}
	}
	fclose(fp_fs);

}

void write_band(complex_3D self)
{
	FILE *fp6;
	int *band_0;
	int ix, iy;
	double kx, ky, ky0, min_band;

	band_0 = ivector(0, nmom);

	fp6 = fopen("mat6.txt", "w");



	for (ix = 0; ix <= nmom; ix++)
	{
		kx = step*(ix);
		band_0[ix] = 0;
		min_band = 430;

		for (iy = 0; iy <= nmom; iy++)
		{
			ky = step*(iy);

			if (min_band >= fabs(xi_band(kx, ky) + self.x1[0][ix][iy]))
			{
				min_band = fabs(xi_band(kx, ky) + self.x1[0][ix][iy]);
				band_0[ix] = iy;
			}
		}
		ky0 = step*(band_0[ix]);
		fprintf(fp6, "%d\t%d\t%f\t%f\t%f\n", ix, band_0[ix], kx, ky0, 180 * atan((1.0 - ky0) / (1.0 - kx)) / M_PI);
	}
	fclose(fp6);

}


/*


void write_band_self(complex_3D self)
{
	FILE *fp31;
	char fi_se4[40], fid4[20];
	int *band_0;
	int ix, iy, ie, w;
	double kx, ky, ky0, min_band, ener;

	band_0 = ivector(0, nmom);

	sprintf(fid4, "dat3");

	mkdir(fid4, 1017);


	for (ix = 0; ix <= nmom; ix++)
	{
		kx = step*(ix);
		band_0[ix] = 0;
		min_band = 430;

		for (iy = 0; iy <= nmom; iy++)
		{
			ky = step*(iy);

			if (min_band >= fabs(xi_band(kx, ky) + self.x1[0][ix][iy]))
			{
				min_band = fabs(xi_band(kx, ky) + self.x1[0][ix][iy]);
				band_0[ix] = iy;
			}
		}
		iy = band_0[ix];
		ky = step*(iy);
		
		sprintf(fi_se4, "dat3/%d.dat", ix);
		fp31 = fopen(fi_se4, "w");
		for (ie = 1; ie <= Nener; ie++)
		{
			ener = energy_fermion(ie);
			fprintf(fp31, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ener, self.s1[ie][ix][iy], self.s2[ie][ix][iy], self.x1[ie][ix][iy], self.x2[ie][ix][iy], self.p1[ie][ix][iy], self.p2[ie][ix][iy], kx, ky, 180 * atan((1.0 - ky) / (1.0 - kx)) / M_PI);

		}


	}
	fclose(fp31);

}
*/


void write_band_self(complex_3D self, double_3D A)
{
	FILE *fp31;
	char fi_se4[40], fid4[20];
	int *band_0;
	int ix, iy, ie, w;
	double kx, ky, ky0, min_band, ener;

	band_0 = ivector(0, nmom);

	sprintf(fid4, "dat3");
	
	mkdir(fid4, 1017);


	for (ix = 0; ix <= nmom; ix++)
	{
		kx = step*(ix);
		band_0[ix] = 0;
		min_band = 430;

		for (iy = 0; iy <= nmom; iy++)
		{
			ky = step*(iy);

			if (min_band >= fabs(xi_band(kx, ky) + self.x1[0][ix][iy]))
			{
				min_band = fabs(xi_band(kx, ky) + self.x1[0][ix][iy]);
				band_0[ix] = iy;
			}
		}
		iy = band_0[ix];
		ky = step*(iy);
	
		sprintf(fi_se4, "dat3/%d.dat", ix);
		fp31 = fopen(fi_se4, "w");
		for (ie = 1; ie <= Nener; ie++)
		{
			ener = energy_fermion(ie);
			fprintf(fp31, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ener, self.s1[ie][ix][iy], self.s2[ie][ix][iy], self.x1[ie][ix][iy], self.x2[ie][ix][iy], self.p1[ie][ix][iy], -self.p2[ie][ix][iy], kx, ky, 180 * atan((1.0 - ky) / (1.0 - kx)) / M_PI);

		}


	}
	fclose(fp31);

}
	
	









#endif
