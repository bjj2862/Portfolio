#ifndef _para_H
#define _para_H 


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
double delt, deltkk;
double wc = 3 ;
int iwbs = 1, iwbc = 400  ;
int N_iter, N_iter2 ;

double mixing_rate_eff, mixing_rate_imp ;

void set_parameters()
{

	enstep = (1.0*wc)/Nener ;
	small = 3.0*enstep/1.0 ;
	step = 1.0/nmom ;
	step2 = 1.0/nmom5 ;
	delt = -3.0*enstep/1.0 ;
	deltkk = enstep/10000.0 ;

	mixing_rate_eff = 0.5;
	mixing_rate_imp = 0.5;
	N_iter = 15;
	N_iter2 = 0;

}



#endif
