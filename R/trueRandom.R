# truerandom.R
#
# R methods to obtain true random numbers from the /dev/random
# entropy collector.
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

#
# Uses Linux /dev/random if available to get high-quality true random
# numbers for use in setting the PRNG seed. Uses hotbits to get
# numbers as well.
#


"runifS" <- function(n, ... , period=10000) {
	k = floor(n/period)
	resetSeed()
	res=runif(n-(period*k), ...)
	if (k>0) {
	   for (i in 1:k) {
		resetSeed()
		res = c(res,runif(period))
	   }
	}
	return(res)
}

"runifT" <- function(n, min=0, max=1) {
	if (min>=max) {
		stop("Max must be > min")
	}
	tmp = trueRandom(n)
	if (is.null(tmp)) {
		warning("No entropy available, returning pseudo-random numbers")
		r = runif(n, min, max)
	} else { 
		r =(tmp/.Machine$integer.max  + 1) * ((max-min)/2)
	}
	return(r)
}

"trueRandom" <-function (n) {
	if (length(n)>1) {
		size=length(n)
	} else {
		size = n
	}

	if (!exists(".EntropyPool",envir=.GlobalEnv)) {
		pool=refreshPool(silent=TRUE)
		if (is.null(pool)) {
			return(NULL)
		}
	} else {
		pool=get(".EntropyPool",envir=.GlobalEnv)
	}

	tr = integer(size)
	i = 1
	while (i<=size) {
		if (pool$current < 1)  {
			pool=refreshPool(silent=TRUE)
			next
		} 
		tr[i] = pool$pool[pool$current]
		pool$current = pool$current -1
		i= i+1
	}	

	assign(".EntropyPool",pool,envir=.GlobalEnv)
	return(tr)
}

"resetSeed" <-function() {
	s = trueRandom(1)
	if (is.null(s)) {
		warning("No entropy available, using system time as seed.")
		s = as.integer(Sys.time())
	}
	set.seed(s)
}

"initPool"<-function(size=512, hbok=TRUE, devrndok=TRUE, silent=FALSE) {
	entropypool = list()
	entropypool$size=size
	entropypool$current=0
	entropypool$pool=integer(length=size)
	class(entropypool)="EntropyPool"
	
	if (size<1 || size > 10000/.Machine$sizeof.long) {
		warning("size out of range")
		return (NULL)
	}

	w = options("warn")
	options(warn=-1)
	if (devrndok) {
	   tri = integer()
           tr = try({tri=readBin('/dev/random', integer(0), signed=FALSE)},
                silent=TRUE)
           if (inherits(tr, "try-error") || (length(tri) == 0)) {
                devrndok=FALSE
           }   
	}
  
	if (hbok) {
	   tri = integer()
	   hb = try({hburl(bytes= .Machine$sizeof.long)})
           if (is.null(hb) || inherits(hb, "try-error")) {
                hbok=FALSE
	   } else {
           	tr = try({tri=readBin(hb, integer(0), signed=FALSE)}, silent=TRUE)
                 if (inherits(tr, "try-error") || (length(tri) == 0)) {
                   hbok=FALSE
                 }  else {
	   		try(close(hb), silent=TRUE)
		}
	   }
	}
	options(warn=as.integer(w))
	entropypool$hbok=hbok
	entropypool$devrndok=devrndok
	if (!devrndok && !hbok ) {
		if (!silent) {
			warning("initialization failed, no true random sources found")
		}
		return(NULL)
	}
	assign(".EntropyPool",entropypool,envir=.GlobalEnv)
	return(entropypool)
}


"hburl" <-function(bytes=1,fmt="bin") {
	hbstring = paste (
		"http://www.fourmilab.ch/cgi-bin/uncgi/Hotbits?"
		,"nbytes=",bytes,"&fmt=",fmt, sep="")
	return(url(hbstring,open="rb"))
}

"refreshPool"<-function(silent=FALSE) {
	if (!exists(".EntropyPool",envir=.GlobalEnv)) {
		if (is.null(initPool(silent=TRUE))) { 
			if (!silent) {
				warning("Could not create pool")
			}
			return(NULL)
		}
	}
	pool=get(".EntropyPool",envir=.GlobalEnv)
	if (pool$devrndok) {
		con = "/dev/random"
		pool$pool = readBin(con, integer(0), signed=FALSE,n=pool$size)
	} else if (pool$hbok) {
		con = hburl(bytes=pool$size)
		pool$pool = readBin(con, integer(0), signed=FALSE,n=pool$size)
		close(con)
	}
	pool$current=length(pool$pool)
	assign(".EntropyPool",pool,envir=.GlobalEnv)
	return(pool)
}
