/*
*  frexp.c
*
* R wrapper for frexp()
*
* Part of the Accuracy package. Available from www.r-project.org and
* www.hmdc.harvard.edu/numerical_issues/
*
*    Copyright (C) 2004  Micah Altman
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

#include <math.h>
#include <time.h>

void R_frexp(double v[], int *nv, double mantissa[], int exp[]) {

	int i, tmpi;
	                                                         
	for (i=0; i < *nv; i++) {
		mantissa[i] = 	frexp(v[i], &tmpi);
		exp[i]=tmpi;
	}
}

void R_time (int *seconds) {
    *seconds = (int) time(NULL);
}
  