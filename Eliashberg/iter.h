#ifndef _iter_H
#define _iter_H


double iter_eff_self(FILE *fp, complex_3D self_eff, double_3D A, double ***Vq_eff)
{
	int ie, je, ix, iy, iter1 ;
	double *difference, difference_0 ; 
	complex_3D self_eff_new;

	self_eff_new = initialize_complex();

// fourier den start
	for(ie=1;ie<=Nener;ie++) 
	{
		fourn_m_to_r(A.s[ie],1);
		fourn_m_to_r(A.x[ie],1);
		fourn_m_to_r(A.p[ie],-1);
	}
// fourier den end
	cal_eff_self_energy(self_eff_new,A,Vq_eff);

	for(ie=1;ie<=Nener;ie++) 
	{
		fourn_r_to_m(self_eff_new.s1[ie],self_eff_new.s2[ie],1);
		fourn_r_to_m(self_eff_new.x1[ie],self_eff_new.x2[ie],1);
		fourn_r_to_m(self_eff_new.p1[ie],self_eff_new.p2[ie],-1);
	}

	difference = cal_difference(self_eff, self_eff_new);
	difference_0 = difference[0];
	printf("difference : %.8f\t%.8f\t%.8f\t%.8f\n", difference[0],difference[1],difference[2],difference[3]) ;
	fprintf(fp, "difference : %.8f\t%.8f\t%.8f\t%.8f\n", difference[0],difference[1],difference[2],difference[3]) ;

	free_dvector(difference,0,3);
	mixing(self_eff, self_eff_new,mixing_rate_eff);
	free_complex_3D(self_eff_new);

	return difference_0;
}



double iter_imp_self(FILE *fp, complex_3D self_imp, complex_3D G, double **Vq_imp)
{
	int ie, je, ix, iy ;
	double *difference, difference_0 ; 
	complex_3D self_imp_new;

	self_imp_new = initialize_complex();

// fourier den start
	for(ie=1;ie<=Nener;ie++) 
	{
		fourn_m_to_r2(G.s1[ie],G.s2[ie],1);
		fourn_m_to_r2(G.x1[ie],G.x2[ie],1);
		fourn_m_to_r2(G.p1[ie],G.p2[ie],-1);
	}
// fourier den end
	cal_imp_self_energy(self_imp_new,G,Vq_imp);

	for(ie=1;ie<=Nener;ie++) 
	{
		fourn_r_to_m(self_imp_new.s1[ie],self_imp_new.s2[ie],1);
		fourn_r_to_m(self_imp_new.x1[ie],self_imp_new.x2[ie],1);
		fourn_r_to_m(self_imp_new.p1[ie],self_imp_new.p2[ie],-1);
	}

	difference = cal_difference(self_imp, self_imp_new);
	difference_0 = difference[0];
	printf("difference : %.8f\t%.8f\t%.8f\t%.8f\n", difference[0],difference[1],difference[2],difference[3]) ;
	fprintf(fp, "difference : %.8f\t%.8f\t%.8f\t%.8f\n", difference[0],difference[1],difference[2],difference[3]) ;

	free_dvector(difference,0,3);
	mixing(self_imp, self_imp_new,mixing_rate_imp);
	free_complex_3D(self_imp_new);

	return difference_0;
}

#endif
