/*
*  randread.c
*
* Read /dev/random with proper timeouts since R nonblocking i/o 
* is not reliable
*
*
*
* Part of the Accuracy package. Available from www.r-project.org and
* www.hmdc.harvard.edu/numerical_issues/
*
*    Copyright (C) 2006  Micah Altman
*
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program; if not, write to the Free Software
*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <fcntl.h>
#ifdef MSVC
#include <io.h>
#define open _open
#define read _read
#define close _close
#else
#include <unistd.h>
#endif

#ifdef NONONBLOCK
#ifdef HASALARM
#include <setjmp.h>
#include <signal.h>
 jmp_buf jumpenv;
#endif
#endif

#ifdef DEBUG
#include <stdio.h>
#endif

int readrand(unsigned int numWanted,unsigned int maxTries, unsigned int sleepInt, int *buf);

void R_readrand(int *numWanted, int *maxTries, int *sleepInt, int *numRead, int result[]) {
	*numRead = readrand(*numWanted, *maxTries, *sleepInt, result);
	#ifdef DEBUG
	int i;
	for (i=0; i<*numRead; i++) {
	  printf("%d\n",(int)result[i]);
	}
	#endif
}

void R_intinfo(int *sizei) {
	*sizei = sizeof(int);
}

#ifdef NONONBLOCK
#ifdef HASALARM
void catch_alarm(int sig_num)
{
         longjmp(jumpenv, 1);
}
#endif
#endif

int readrand (unsigned int numWanted ,unsigned int maxTries,unsigned int sleepInt,
	int *result) {
   
   int rv=0, bytesRead=0,  bytesWanted = numWanted * sizeof(int);
   unsigned int numTries=0; 
   #ifdef NONONBLOCK
   int fd=open("/dev/random",O_RDONLY);
   #ifdef HASALARM
   signal(SIGALRM, catch_alarm);
   alarm(sleepInt*maxTries);
   if (setjmp(jumpenv) != 0) {
        return(-1);
   }
   #endif
   #else
   int fd=open("/dev/random",O_NONBLOCK|O_RDONLY);
   #endif
   
  
   if (fd<0) {
    	return (-1);
   }
   

  
   while (bytesRead<bytesWanted && numTries<maxTries) {
   	rv = read(fd, (((unsigned char*) result)+bytesRead), bytesWanted-bytesRead);
	#ifdef DEBUG
	printf ("read %d \n", rv);
	#endif
	if (rv>0) {
		bytesRead+=rv;
	}
	#ifndef NONONBLOCK
	#ifdef HASALARM
	if (bytesRead<bytesWanted && numTries<maxTries) {
  		sleep(sleepInt);	 
	}
	#endif
	#endif
	numTries++;
   }
   #ifdef NONONBLOCK
   #ifdef HASALARM
   alarm(0);
   #endif
   #endif
   close(fd);
   return(bytesRead/sizeof(int));
}

#ifdef DEBUG
main() {
	int intsWanted=20, maxTries=2, sleepInt=2, intsRead=0; 
	int result[20];	

	
	R_readrand(&intsWanted, &maxTries, &sleepInt, &intsRead, result);
	
	printf("Read %d \n", intsRead);
	int i;
	for (i=0; i < intsRead; i++) {
		printf (":%d\n",result[i]);
	}
}
#endif
