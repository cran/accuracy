######################################################
#       
# accFirstLook
#       
# Provides an overview first look at a dataset.
#
# For exploratory analysis.
#
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################


firstLook<-function(x,forceType=TRUE) {
	ret <- list()
	if (!is.data.frame(x) && !is.ts(x)) {
		x<-as.data.frame(x)
	}
	if (forceType && !is.ts(x)) {
		notnum<-which(!sapply(x,is.numeric))
		x[notnum] <-	sapply(x[notnum], simplify=FALSE,as.factor)
	}
	ret$data<-x
	
	nvar<-dim(x)[2]
	if (is.null(nvar)) nvar<-1

	if (nvar>1) {
	  # pairwise.complete.obs , with Kendall would be better
    # but is very slow with kendall method
		ret$cor <- cor(x,use="complete.obs", method="spearman")
		pc<-NULL
    if (is.data.frame(x)) {
      numeric.columns <- which(sapply(x,is.numeric))
	  	if (length(numeric.columns)>1)  {
			  pc <- try(princomp(x[numeric.columns], cor=TRUE),silent=T)
	  	}
 	   }
 	   if (is.ts(x)&&(nvar>1)) {
 	      pc <-try(princomp(x,cor=TRUE),silent=T)
     }
     if(!is.null(pc) && !inherits(pc,"try-error")){
      ret$pc<-pc
     }
	}
	class(ret)<-"accFirstLook"
	return (ret)
}


print.accFirstLook<-function(x,...) {
  nvar<-dim(x$data)[2]
	if (is.null(nvar)) nvar<-1
   if (nvar==1) {
		if (is.numeric(x$data[[1]])) {		
			print(summary(x$data[[1]]),...)
		} else {
		  	tmp.t <- table(x$data[[1]])
			print(tmp.t,...)
			print(summary(tmp.t,...))
		}
	} else {
	  if (is.data.frame(x$data)){
	    isfact<-sapply(x$data,is.factor)
      x$data[isfact] <-	sapply(x$data[isfact],as.numeric)
    }
	  print(summary(x$data)) 
		if (!is.null(x$pc)) {
			print(summary(x$pc),...)
		}
		print(x$cor,...)
	}
}
	
plot.accFirstLook<-function(x,single=FALSE, ask=dev.interactive(),...) {
   op <- par(no.readonly=T)
   on.exit(par(op))
   par(ask = ask)
 
  nvar <- dim(x$data)[2]
 	if (is.null(nvar)) nvar<-1
	if (nvar==1) {
	  plot.accSingle(x$data)
	} else {
	  heatmap(x$cor,margins=c(10,10),revC=T,scale="none",...)
		if (!single) {
		   if (is.ts(x$data)) {
          plot(x$data,main="")
		   }
		   plotvar<-min(nvar,25)
		   row2<-ifelse(plotvar>floor(sqrt(plotvar))*ceiling(sqrt(plotvar)),ceiling(sqrt(plotvar)),
		      floor(sqrt(plotvar)))
       par( mfrow=c(ceiling(sqrt(plotvar)),row2))
       for (i in 1:plotvar) {
          if (is.ts(x$data)) {
             plot.accSingle(as.data.frame(x$data[,i]))
          } else {
             plot.accSingle(x$data[i])
          }
       }
       par( mfrow=c(1,1)) 
       if (is.data.frame(x$data)) {
		     isfact<-sapply(x$data,is.factor)
	 	     x$data[isfact] <-	sapply(x$data[isfact],as.numeric)
		     numeric.data <- x$data[which(sapply(x$data,is.numeric))]

       } else {
         numeric.data<-x$data
       }
       nvar.num<-dim(numeric.data)[2]
       if (nvar.num>1) {
			     ow<-options(warn=0)
			     .require<-require                
			     if (.require(quietly=T,"MASS"))  {
			       #par(mfrow=c(nvar.num,2))
			       #for (i in 1:nvar.num) {
             #   i=1
				     #   parcoord(numeric.data[c((1:nvar.num+i)%%nvar.num+1)])
		         #}
		          parcoord(numeric.data)
			     }
			       par( mfrow=c(1,1)) 
			     options(ow)
       }                                           
 

		   if (!is.null(x$pc)) {
			    biplot(x$pc,cex=c(.4,2),...)
		   }
	  }
	}
}

plot.accSingle<-function(x,...) {
    if (is.ts(x)) {
      plot(x)
    } else  if (is.numeric(x[[1]])) {		
			hist(x[[1]],xlab="",main=names(x),...)
		} else {
			mosaicplot(table(x[[1]]),main=names(x),...)
	  }
}
