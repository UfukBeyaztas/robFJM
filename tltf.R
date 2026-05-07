tltf <- function(long, surv, eps_time = 1e-10) {
  
  surv <- surv[order(surv$ID), , drop = FALSE]
  
  t_obs <- surv$time
  names(t_obs) <- as.character(surv$ID)
  
  keep <- long$obstime <= (t_obs[as.character(long$ID)] + eps_time)
  keep[is.na(keep)] <- FALSE
  long2 <- long[keep, , drop = FALSE]
  
  ids <- sort(intersect(unique(long2$ID), unique(surv$ID)))
  long2 <- long2[long2$ID %in% ids, , drop = FALSE]
  surv2 <- surv[surv$ID %in% ids, , drop = FALSE]
  
  long2 <- long2[order(long2$ID, long2$obstime), , drop = FALSE]
  surv2 <- surv2[order(surv2$ID), , drop = FALSE]
  
  list(long = long2, surv = surv2)
}