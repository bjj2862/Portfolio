
#ifndef _STRUCT_H
#define _STRUCT_H

typedef struct
{
	double ***s1, ***s2, ***x1, ***x2, ***p1, ***p2 ; 
} complex_3D;

typedef struct
{
	double *s1, *s2, *x1, *x2, *p1, *p2 ; 
} complex_1D;

typedef struct
{
	double ***s, ***x, ***p;
} double_3D;


typedef struct 
{
	double *delta, *kapa, *xamp;
} Vignolle;



#endif 
