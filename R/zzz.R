.onLoad <- function(lib, pkg) {
  if((version$major<1) ||(version$minor<6 && version$major==1)) {
    stop("This version for R 1.6 or later")
  }
  if (!frexpTest(FALSE)) {
  	warning("frexp: failed self-test");
  }
  if (!perturbTest()) {
  	warning("perturb: failed self-test");
  } 
  if (!secholTest()) {
  	warning("sechol: failed self-test");
  }
  return(TRUE);
}
