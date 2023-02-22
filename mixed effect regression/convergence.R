mod_check <- function(mod) {
  if(is.null(unlist(mod@optinfo$conv$lme4))) {
    retval = 0 #returns 0 no convergence issues (convergence list empty)
  }
  else {
    if (isSingular(mod)) {
      retval = 1 #returns 1 if the model was singular
    } else {
      retval = 2 #returns 2 if there was another convergence problem
    }
  }
  retval #returns the value (1, 0, or -1)
}
