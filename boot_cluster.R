library(boot)
library(plyr)

boot.cluster <- function(data, cluster_var, model, extractor = identity, strata_var = NULL, sub_sample_n = NULL, ...){
  # Add rows to ensure that all clusters have same number of observations
  # Could be optimized, is however only run once
  data[[cluster_var]] <- as.character(data[[cluster_var]])
  while(TRUE){
    observations_per_cluster <- table(data[[cluster_var]])
    if(all(max(observations_per_cluster) == observations_per_cluster)){
      break
    }
    missing_row <- data.frame(names(observations_per_cluster[observations_per_cluster < max(observations_per_cluster)][1]))
    names(missing_row) <- cluster_var
    data <- rbind.fill(data, missing_row)
  }

  orig_data <- data
  
  # Check that all clusters have equal number of observations
  if(!all(observations_per_cluster[1] == observations_per_cluster)){
    stop("All clusters must have same number of observations")
  }
  # Get unique strata values
  stratas <- if(!is.null(strata_var)){unique(orig_data[[strata_var]])} else {NULL}
  
  # Get unique cluster names
  all_clusters <- data.frame(unique(orig_data[c(cluster_var, strata_var)]))
  
  strata_count <- ifelse(is.null(strata_var), 1, length(stratas))
  sample_size <- ifelse(is.null(sub_sample_n), nrow(all_clusters), sub_sample_n * strata_count)
  id_col <- sort(rep(1:sample_size, observations_per_cluster[1]))
  
  sampler <- function(data, index){
    clusters <- cbind(data[index, , drop = FALSE], row = 1:length(index))
    if(!is.null(sub_sample_n)){
      if(!is.null(strata_var)){
        clusters <- ddply(clusters, strata_var, function(x){
          x[1:sub_sample_n, ]
        })
      }
      else{
        clusters <- clusters[1:sub_sample_n, ]
      }
    }
    
    sample_long <- join(clusters, orig_data, by = cluster_var)
    sample_long[["cluster_id"]] <- sample_long[[cluster_var]]
    sample_long[[cluster_var]] <- id_col
    sample_long
  }
  if(!is.null(strata_var)){
    strata <- all_clusters[[strata_var]]
  }
  else{
    strata <- rep(1, nrow(all_clusters))
  }
  boot.custom(data = all_clusters, sampler = sampler, model = model, extractor = extractor, strata = strata, ...)
}

boot.custom <- function(data, sampler, model, extractor, ...){
  modelX <- error.safe(model, verbose_error = T)
  statistic <- function(data, index){
    x <- sampler(data, index)
    m <- modelX(x)
    if(is.null(m)){
      return(NA)
    }
    extractor(m)
  }
  boot(data = data, statistic = statistic, ...)
}

error.safe <- function(fun, print_error = TRUE, verbose_error = FALSE, on_fail = NULL){
  function(...){
    res <- try(suppressWarnings(fun(...)), silent = TRUE)
    if (is.element("try-error", class(res))) {
      if(print_error){
        if(verbose_error){
          print(res)
        }
        else{
          print("error")
        }
      }
      on_fail
    } else {
      res
    }
  }
}

count.calls <- function(fun){
  count <- 0
  function(...){
    count <<- count + 1
    print(count)
    fun(...)
  }
}

format2 <- function(x){
  format(x, digits = 2, nsmall = 2)
}

summary.boot.power <- function(bootout, ...){
  p_index <- which(names(bootout$t0) == "p")
  cat.ln("Power = ", format(sum((bootout$t[, p_index] <=.05)) / bootout$R * 100, digits = 4, nsmall = 2), "%")
  for(index in setdiff(1:length(bootout$t0), grep("^p$", names(bootout$t0)))){
    vc <- bootout$t[, index]
    v <- names(bootout$t0)[index]
    cat.ln(v, " mean = ", format2(mean(vc)))
    cat.ln(v, " sd = ", format2(sd(vc)))
    cat.ln(v, " median = ", format2(median(vc)))
    cat.ln(v, " CI = [", format2(quantile(vc, .025)), ", ", format2(quantile(vc, .975)), ']')
  }
  
  dots <- list(...)
  if(length(dots) > 0){
    boot.ci(bootout, ...)
  }
}


boot.bcaCIs <- function(boots){
  res <- list()
  for(i in 1:length(boots$t0)){
    if(!is.na(boots$t0[i])){
      ci.res <- boot.ci(boots, index = i, type = c("bca"))
      x <- tail(ci.res[["bca"]][1, ], 2)
      low = x[1]
      high = x[2]
      res[[names(boots$t0)[i]]] <- c(low, high)
    }
  }
  res
}


boot.CIs <- function(boots, CI){
  types = c("norm" = "normal",
            "basic" = "basic",
            "stud" = "student",
            "perc"= "percent",
            "bca" = "bca")
  res <- list()
  for(i in 1:length(boots$t0)){
    ci.res <- boot.ci(boots, index = i, type = CI)
    x <- tail(ci.res[[types[CI]]][1, ], 2)
    low = x[1]
    high = x[2]
    res[[names(boots$t0)[i]]] <- c(low, high)
  }
  res
}


# Randomness testing
# 
# ids <- NULL
# collect_ids <- function(data){
#   ids <<- c(ids, data$cluster_id)
#   nrow(data)
# }
# 
# boot_res2 <- boot.cluster(malin, "StudyID", collect_ids, count.calls(identity), R = 5000)
# x <- table(ids)
# sum(x) / nrow(malin)
# plot(density(x))
# 
# x1 <- table(ids1)
# x2 <- table(ids2)
# x3 <- table(ids3)
# 
# cor.test(x1[order(names(x1))], x2[order(names(x2))])
# cor.test(x1[order(names(x1))], x3[order(names(x3))])
# cor.test(x2[order(names(x2))], x3[order(names(x3))])