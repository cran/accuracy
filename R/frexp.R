# frexp
#
# Extract mantissa's and exponents from a vector of numbers
# R wrapper around the frexp() C-library call.
#
# Part of the Accuracy package. Available from www.r-project.org and
# www.hmdc.harvard.edu/numerical_issues/
#
#    Copyright (C) 2004  Micah Altman
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


"frexp" <-
function(v) {
	r = replicate(length(v), NA);
	r = cbind(r,r);
	i = which(is.finite(v) | !is.na(v));
	tmp = .C("R_frexp", PACKAGE="accuracy", as.double(v[i]), as.integer(length(v[i])),
		mantissa=double(length(v[i])),
		exponent=integer(length(v[i])));
	r[i,] = cbind(tmp$mantissa,tmp$exponent);
	dimnames(r)[[2]]=c("Mantissa","Exponent")
	return(r);
}


"frexpTest" <-
function(silent=TRUE) {
	d=options()$digits;
	options(digits=15);
	f=round(frexp(c(1,2,4,1.1,2.1,4.1,1000,20000,1.02E213)),digits=5);
	x=cbind(rbind(.5,.5,.5,.55,.525,.5125,.97656,.61035,.75747),
		 rbind(1,2,3,1,2,3,10,15,708));
	options(digits=d);
	ret=(sum(f==x)==18);
	if (!ret && !silent) {
		warning("Failed frexp self test.");
	}
	return(ret);
}

"LRE"<-
function(x,correct, use.LAE=TRUE) {
        zeros = which(x==0);
        nonzeros = which(x!=0);
        nas = which(x==NA);

        res=double(length(x));

        res[nas]=NA;
        res[nonzeros]= -1 * 
                log10(abs(x[nonzeros]-correct[nonzeros])/correct[nonzeros]); 
        if (use.LAE) {
                res[zeros]=-1*log10(abs(x[zeros]-correct[zeros]));
        } else {
                res[zeros]=NA;
        }

        return(res);
}

