HTML.perturbS<-function(x, ...) {
   HTML("Sensitivity of coefficients to perturbations:\n",...)
   HTML(hSummary(attr(x,"coef.betas.m")),...)

    if (!is.null(attr(x,"coef.stderrs.m"))) {
        HTML ("\n\nSensitivity of stderrs to perturbations:\n",...)
        HTML (hSummary(attr(x,"coef.stderrs.m")),...)
     }
}

HTML.perturbAnova<-function(x,...) {
        HTML("\nsumsq:\n\n",...)
        HTML(hSummary(attr(x,"coef.sumsq.m")), ...)
        HTML("\nmeansq:\n\n",...)
        HTML(hSummary(attr(x,"coef.meansq.m")),...)
        HTML("\npr:\n\n",...)
        HTML(hSummary(attr(x,"coef.pr.m")),...)

}
