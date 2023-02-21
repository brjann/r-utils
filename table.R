table.decimals <- function(v, d=2){
  format(v, nsmall = d, digits = d, scientific = F)
}

table.format3p <- function(v, add.equals = T){
  if (round(v, 3) == 0){
    return ("< .001")
  }
  return (paste0(ifelse(add.equals, "= ", ""), substring(table.decimals(v, 3), 2)))
}

table.pstars <- function(p, p.10 = FALSE){
  if(p <= .001){
    '***'
  }
  else if(p <= .01){
    '**'
  }
  else if(p <= .05){
    '*'
  }
  else if(p.10 & p <= .10){
    "."
  }
}

table.mean.se <- function(mean, se){
  paste0(table.decimals(mean), " (", table.decimals(se), ")")
}

table.CI <- function(CIs){
  CI_l <- min(CIs)
  CI_u <- max(CIs)
  if(!is.null(CIs)){
    paste0("[", table.decimals(CI_l), ", ", table.decimals(CI_u), ']')
  }
  else{
    "[xx, xx]"
  }
}
