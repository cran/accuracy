# truerandom.R
#
# R methods to obtain true random numbers from the /dev/random
# entropy collector and web based sources
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
# These routines use the Linux /dev/random if available to get high-quality true random
# numbers for use in setting the PRNG seed. Uses hotbits to get
# numbers as well.
#
# For efficiency, the entropy gathered from these sources is cached in a local 'pool' of entropy,
# and only refreshed when the pool is empty. Each of the exposed routines draws from
# the same pool


######################################################
#
# rstarMix
#
# Calls random number generator reseeding periodically from entropy sources
# 
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################


"rstarMix" <- function(n, ... , period=10000,maxTries=10,silent=TRUE,rfunc="runif") {
	res<-numeric(n)
	k = floor(n/period)
	resetSeed(maxTries=maxTries,silent=silent)
	remainder <- n-(period*k)
	if (remainder>0) {
		res[1:remainder] <- do.call(rfunc, list(n=remainder,...))
	}
	if (k>0) {
	   for (i in 1:k) {
		resetSeed(maxTries=maxTries,silent=silent)
		res[(remainder+((i-1)*period)+1):(remainder+(i*period))] <-
			 do.call(rfunc, list(n=period,...))
	   }
	}
	return(res)
}

######################################################
#
# runifS
#
# Calls runif() reseeding periodically from entropy sources
# 
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################


"runifS" <- function(n, ... , period=10000,maxTries=10,silent=TRUE) {
	rstarMix(n=n,...,period=period,maxTries=maxTries,silent=silent)
}
"rnormS" <- function(n, ... , period=10000,maxTries=10,silent=TRUE) {
	rstarMix(n=n,...,period=period,maxTries=maxTries,silent=silent,rfunc="rnorm")
}

######################################################
#
# runifT
#       
# Returns random uniform values drawn from entropy pools
# 
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################


"runifT" <- function(n, min=0, max=1,maxTries=10,silent=TRUE) {
	if (min>=max) {
		stop("Max must be > min")
	}
	tmp <- trueRandom(n,maxTries=maxTries,silent=silent)
	#r <--((((tmp/.Machine$integer.max)+1)/2)*(max-min)) + min
	r <-(as.numeric(tmp)+.Machine$integer.max)/(2*.Machine$integer.max) * (max-min) + min
	if (length(tmp)<n) {
		if(!silent) {	
			warning("Not enough entropy available, returning some pseudo-random numbers")
		}
		r = c(r,runif(n-length(r), min, max))
		sample(r,n)
	}
	return(r)
}

######################################################
#
# trueRandom
#       
# Returns raw random bits
# 
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
#
######################################################


"trueRandom" <-function (n, maxTries=100,silent=TRUE) {
	if (length(n)>1) {
		size=length(n)
	} else {
		size = n
	}

	tr=integer(0)
  pool=getPool(silent=silent)

	tries=0;
	while ( (tries<maxTries) && (length(tr)<size)) {
		if (pool$current < 1)  {
			pool=refreshPool(silent=silent)
			tries = tries+1
			if ((tries<maxTries) && (pool$current<1)) {
			  
			  if (is.null(options()$timeout)) {
			      timeout = 60
        }  else {
            timeout = options()$timeout
        }
				Sys.sleep(timeout)
			}
			next
		} 
		needed = size-length(tr)
		chunk = min(needed,pool$current)
		tmp = pool$pool[(pool$current-chunk+1):pool$current]	
		tr = c(tr,tmp)
		pool$current = pool$current-chunk
		setPool(pool)
	}
	return(tr)
}

######################################################
#
# resetSeed
#       
# Resets seed with true random value
# 
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################


"resetSeed" <-function(maxTries=10,silent=TRUE,kind=NULL,normal.kind=NULL) {
	oldSeed = NULL
	# deal with random number kind setting
	if (!is.null(kind)) {
		RNGkind(kind=kind)
	}
	if (!is.null(normal.kind)) {
		RNGkind(kind=normal.kind)
	}
	
	# generate a number so that .Random.seed gets set
	if (exists(".Random.seed")) {
	   oldSeed=get(".Random.seed")
  	}
	runif(1)	
  	s = trueRandom(length(.Random.seed),maxTries=maxTries,silent=silent)
	
	if (is.null(s)|| length(s)==0) {
		if (!silent) {
			warning("No entropy available, using system time as seed.")
		}
    s=Sys.time() 
  }
	retval<-set.seed(s,kind=kind,normal.kind=normal.kind)
	invisible(retval)
}

######################################################
#
# initPool
#       
# Initializes the entropy pool storage structure
# 
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################



"initPool"<-function(hbok=0, rook=0, devrndok=0, silent=FALSE) {
	entropypool = list()
	entropypool$current=0
	entropypool$pool=integer(0)
	entropypool$intinfo = 
	  .C("R_intinfo", PACKAGE="accuracy", sizeof.int =integer(1),  DUP=FALSE)
	class(entropypool)="EntropyPool"



   if (!file.exists("/dev/random")) {
          devrndok=.Machine$integer.max
   } else {
           devrndok=0
   }
	entropypool$hbok=hbok
	entropypool$rook=rook
	entropypool$devrndok=devrndok

	return(entropypool)
}

######################################################
#
# hbUrl
#
# [Internal Function]
#       
# generates the hotbytes url
# 
# Parameters:
#
# - bytes --  number of bytes requested
# - fmt -- hotbits return value (use "bin" except
#	for debugging!)
# 
######################################################


"hburl" <-function(bytes=256,fmt="bin") {
	bytes=round(bytes)
	if (bytes<1) {
		bytes =1
	} else if (bytes>2048) {
		bytes=2048
	}
	hbstring = paste (
		"http://www.fourmilab.ch/cgi-bin/uncgi/Hotbits?"
		,"nbytes=", bytes,"&fmt=",fmt, sep="")
	rtval=tryR(url(hbstring,open="rb"), silent=TRUE)
	if (inherits(rtval,"try-error")) {
		return(NULL)
	} else {	
		return(rtval)
	}
}

"rourl" <-function(bytes=256,fmt="f") {
	bytes=round(bytes)
	if (bytes<1) {
		bytes =1
	} else if (bytes>16384) {
		bytes=16384
	}
	hbstring = paste (
		"http://www.random.org/cgi-bin/randbyte?"
		,"nbytes=", bytes,"&format=",fmt, sep="")
	rtval=tryR(url(hbstring,open="rb"), silent=TRUE)
	if (inherits(rtval,"try-error")) {
		return(NULL)
	} else {	
		return(rtval)
	}
}

######################################################
#
# refreshPool
#
# [Internal Function]
#       
# refreshes the entropy pool storage structure 
# 
# Parameters:
#
# silent - print debugging messages
# 
######################################################

"refreshPool"<-function(silent=FALSE) {

	pool=getPool(silent=silent)
	oldpool=pool$pool
	pool$pool=integer(0)
	st = Sys.time()
	if (is.null(options()$timeout) ) {
	   timeo=60
   } else {
   
	   timeo = options()$timeout
	}

	if (silent) {
		ow= options(warn=-1)
	} else {
	  ow =options("warn")
  }
	if ((st>pool$rook) && (length(pool$pool)==0)) {
	   con = rourl()
	   if (is.null(con)) {
			if (!silent) {warning("RO connection ERROR");}
			pool$rook=st+timeo
	   } else {
           	tmp = tryR({readBin(con, integer(0), signed=FALSE, n=65536)}, silent=TRUE)
	           close(con)
             if (inherits(tmp, "try-error")) {
			if (!silent) {warning("RO ERROR");}
			pool$rook=st+timeo
		} else if (sum(oldpool==tmp)>0) {
			if (!silent) 	{warning("RO REPEATED POOL")}
			pool$rook=st+10*timeo
		} else if (length(tmp)==0) {
			if (!silent)  {warning("RO NO POOL")}
			pool$rook=st+timeo
		} else {
			pool$pool = tmp
		}
	   }
	}

	if ((st>pool$hbok) && (length(pool$pool)==0)) {
	   con = hburl()
	   if (is.null(con)) {
			if (!silent) {warning("HB connection ERROR");}
			pool$hbok=st+timeo
	   } else {
           	tmp = tryR({readBin(con, integer(0), signed=FALSE, n=65536)}, silent=TRUE)
		close(con)

                if (inherits(tmp, "try-error")) {
			if (!silent) {warning("HB COMMUNICATIONS ERROR");}
			pool$hbok=st+timeo
		} else if (isTRUE(all.equal(tmp, 
			c( 544567129, 1702257000, 1668834592, 1701078373,1870209124, 840987253 , 1869098292,1897951861,1635020661,1919903264,1953450016,1937008962 )))) {
			if (!silent)	{warning("HB too many request");}
			pool$hbok=st+180*timeo
		} else if (sum(oldpool==tmp)>0) {
			if (!silent) 	{warning("HB REPEATED POOL")}
			pool$hbok=st+10*timeo
		} else if (length(tmp)==0) {
			if (!silent)  {warning("HB NO POOL")}
			pool$hbok=st+timeo
		} else {
			pool$pool = tmp
		}
	   }
	} 

	if ((st>pool$devrndok) && (length(pool$pool)==0)) {
		tmp= readDevRand()
		if (isTRUE(all.equal(oldpool,tmp))) {
			if (!silent) {
				warning("DEV RANDOM REPEATED POOL")
			}
		}
		pool$pool = tmp
	} 

  options(ow)
	pool$current=length(pool$pool)
	setPool(pool)
	return(pool)
}


######################################################
#
# readDevRand
#
# [Internal Function]
#       
# Safely reads /dev/random
#
# This is a workaround for the fact that /dev/random blocks
# and R does not handle timeouts on blocking I/O correctly
# (in fact, open() can block even when block is set to false :-(
# 
# Parameters:
#
# timeout = timeout value
# numInts = number of ints to read
#
######################################################


readDevRand<-function(numInts=32,timeout=options("timeout")[[1]]) {
	if (!file.exists("/dev/random")) {
		return(integer(0))
	}
	 sleepInt=5
	 maxTries=ceiling(timeout/sleepInt)
	 tmp = .C("R_readrand", PACKAGE="accuracy",  as.integer(numInts),
		as.integer(maxTries), as.integer(sleepInt),
		 numRead =integer(1), results=integer(numInts), DUP=FALSE)
	if (tmp$numRead>0) {
		length(tmp$results)=tmp$numRead
		return(tmp$results)
	} else {
		integer(0)
	}
}


#
#  Internal methods
#

getPool<-function(silent=TRUE) {
  tr=integer(0)

  if (is.R()) {
     existsPool = exists(".EntropyPool",envir=.GlobalEnv)
  } else {
     existsPool = exists(".EntropyPool",where=1)
  }
  
  
	if (!existsPool) {
		pool=initPool(silent=silent)
		if (is.null(pool)) {
			return(tr)
		}
	} else {
	  if (is.R()){
		    pool=get(".EntropyPool",envir=.GlobalEnv)
    } else {
        pool=get(".EntropyPool",where=1)
    }
	}
	return(pool)
}

setPool<-function(pool) {
    if (is.R()){
		    assign(".EntropyPool",pool,envir=.GlobalEnv)
    } else {
        assign(".EntropyPool",pool,where=1)
    }
}
