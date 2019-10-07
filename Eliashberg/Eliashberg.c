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
#include <omp.h>

#include "struct.h"
#include "para.h"

#include "some_function.h"
#include "write.h"

#include "initialize.h"
#include "print_time.h"
#include "cal_function.h"
#include "iter.h"


int main(int argc, char *argv[])
{

//define variable 
	FILE *fp;
	int iter, iter2 ;	
	double ***Vq_eff, **Vq_imp;
	double difference_eff, difference_imp ;
	time_t time_0, time_now, running_time ;

// define struct
	complex_3D self_eff, self_eff_new, self_imp, self_imp_new, G;
	double_3D A;
	complex_1D DOS;


//end define variable

	fp = fopen("output.txt","w") ;

	time(&time_0);
	printf("%s",ctime(&time_0));
	fprintf(fp,"%s",ctime(&time_0));

	set_parameters();





	self_eff = initialize_eff();
	self_imp = initialize_imp();
	G = initialize_complex();
	A = initialize_real();

	Vq_imp = initialize_Vq_imp();
	Vq_eff = initialize_Vq_spin();

	Cal_Green_function(G,A,self_eff,self_imp) ;

// time
	print_time(fp,time_0);

	difference_eff = 1 ;

	for(iter=1;iter<=N_iter;iter++)
	{
		printf("iter1 = %d .\t The difference_eff = %.8f .\n" , iter,difference_eff);
		fprintf(fp,"iter1 = %d .\t The difference_eff = %.8f .\n" , iter,difference_eff);

		difference_eff = iter_eff_self(fp, self_eff,A, Vq_eff);
		Cal_Green_function(G,A,self_eff,self_imp) ;
		print_time(fp,time_0);

		difference_imp = 1 ;

		for(iter2=1;iter2<=N_iter2;iter2++)
		{
			printf("iter2 = %d .\t The difference_imp = %.8f .\n" , iter,difference_imp);
			fprintf(fp,"iter2 = %d .\t The difference_imp = %.8f .\n" , iter,difference_imp);
			difference_imp = iter_imp_self(fp, self_imp,G, Vq_imp);
			Cal_Green_function(G,A,self_eff,self_imp) ;
			print_time(fp,time_0);
			if(difference_imp<=tol)
			{
				iter2 = iter2 + N_iter2;
			}

		}

		if(difference_eff<=tol)
		{
			iter = iter+N_iter;
		}
	}

	DOS = Cal_DOS(G,A);

	write_DOS(DOS);

	write_minus_spectral_function(fp, A);
	write_Fermi_surface(A);
	write_band(self_eff);
	write_band_self(self_eff,A);
//	write_self(self_eff,"self_eff.txt");

	fclose(fp);

}

