#ifndef _print_time_H
#define _print_time_H


#include <time.h>


void print_time(FILE *fp,time_t time_0)
{
	time_t time_now, running_time ;
	time(&time_now);
	printf("%s",ctime(&time_now));
	fprintf(fp,"%s",ctime(&time_now));
	running_time = time_now - time_0   ;
	printf("h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	
	fprintf(fp,"h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);	
}

#endif
