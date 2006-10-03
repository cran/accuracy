# R bridge
#
# Functions to bridge differences in R and Splus 



# simplify uses the SIMPLIFY=TRUE logic from sapply()
# 
# This helps when porting code originally using sapply()
# which won't work in S-plus because of scoping differences


   simplify<-function(answer) {
      if (length(answer) && length(common.len <- unique(unlist(lapply(answer, 
        length)))) == 1) {
        if (common.len == 1) 
            unlist(answer, recursive = FALSE)
        else if (common.len > 1) 
            array(unlist(answer, recursive = FALSE), dim = c(common.len, 
                length(answer)), dimnames = if (!(is.null(n1 <- names(answer[[1]])) & 
                is.null(n2 <- names(answer)))) 
                list(n1, n2))
        else answer
     } 
     else answer
   }

# tryR() provides R style try behavior in R or S-Plus

   tryR<-function(expr,silent=FALSE) {
      if (is.R()){
          return(try(eval.parent( expr),silent))
      } else {
          tmp = try(eval.parent(expr),silent)
          if (inherits(tmp,"Error")) {
              if (!silent) {
                print (tmp)
              }
              class(tmp)="try-error"
          }
          return(tmp)
      }
   }
   
   
# Replacement functions for some R functions not found in SPLUS
#
#   url() and close.url()
#   replicate()
#   isTRUE()
#   Sys.time() -- note not a true replacement, only returns integer time, not Posix struct
#   readBin() -not a true replacement

if (!is.R()){

  dev.interactive<-function(){
    dev.cur()==2
  }
  
  readBin<-function(con, what, n = 1, size = NA, signed = TRUE, endian = .Platform$endian)  {
    if (!is.na(size) ||endian!=.Platform$endian) {  
        warning("Unsupported options")
    }
    ow=options(warn=-1)
    tmp=readRaw(con,what,n)
    options(ow)
    return(tmp)
  }

  Sys.time<-function(){
      .C("R_time",integer(1),COPY=FALSE)[[1]]
  }

    
   replicate<-
   function (n, expr, simplify = TRUE) { 
        res = vector(length=n,mode="list")
        for (i in 1:n) {
          res[[i]] = eval.parent(match.call()$expr)
        }
        if (simplify) {
          res=simplify(res)
        }  else {
          res
        }
   }
   
   isTRUE<-    
   function (x) 
   identical(TRUE, x)

   #
   #
   #
   


   #
   # URL methods
   # 
 
   url<-
   function (description, open = "", blocking = TRUE, encoding = getOption("encoding")) {
     if (is.null(encoding)){
       encoding="native.enc"
     }
     if (open == "") {
       open="r"
     }
     if (open !="r" && open !="rb") {
      warning("Url's can only be open for reading")
     }
    
     tmpfileName <- tempfile("url")
  
     # no explicit error handling, since URL doesn't either
     # rely on download.file to throw an exception
     download.file(description, tmpfileName)
     con = file(tmpfileName, open=open, blocking=blocking)
     class(con) = "url"
     return(con)
  }

   ow=options(warn=-1)
   setClass("url", "file")
   options(ow)


   close.url<-
   function(con,...) {
     class(con)="file"
     status = close(con,...)
     unlink(unclass(con)$description)
     return(status)
  }


   options(warn=-1)
   setMethod("close","url",close.url)
   options(ow)
   
   # favorite test data!
   longley <-
structure(list(GNP.deflator = c(83, 88.5, 88.2, 89.5, 96.2, 98.1, 
99, 100, 101.2, 104.6, 108.4, 110.8, 112.6, 114.2, 115.7, 116.9
), GNP = c(234.289, 259.426, 258.054, 284.599, 328.975, 346.999, 
365.385, 363.112, 397.469, 419.18, 442.769, 444.546, 482.704, 
502.601, 518.173, 554.894), Unemployed = c(235.6, 232.5, 368.2, 
335.1, 209.9, 193.2, 187, 357.8, 290.4, 282.2, 293.6, 468.1, 
381.3, 393.1, 480.6, 400.7), Armed.Forces = c(159, 145.6, 161.6, 
165, 309.9, 359.4, 354.7, 335, 304.8, 285.7, 279.8, 263.7, 255.2, 
251.4, 257.2, 282.7), Population = c(107.608, 108.632, 109.773, 
110.929, 112.075, 113.27, 115.094, 116.219, 117.388, 118.734, 
120.445, 121.95, 123.366, 125.368, 127.852, 130.081), Year = as.integer(c(1947, 
1948, 1949, 1950, 1951, 1952, 1953, 1954, 1955, 1956, 1957, 1958, 
1959, 1960, 1961, 1962)), Employed = c(60.323, 61.122, 60.171, 
61.187, 63.221, 63.639, 64.989, 63.761, 66.019, 67.857, 68.169, 
66.513, 68.655, 69.564, 69.331, 70.551)), .Names = c("GNP.deflator", 
"GNP", "Unemployed", "Armed.Forces", "Population", "Year", "Employed"
), row.names = c("1947", "1948", "1949", "1950", "1951", "1952", 
"1953", "1954", "1955", "1956", "1957", "1958", "1959", "1960", 
"1961", "1962"), class = "data.frame")

   
  
} #END is.R() BLOCK
