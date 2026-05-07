build_designs_std <- function(dat_surv, dat_long, long_formula, surv_formula) {
  dat_surv <- dat_surv[!duplicated(dat_surv$id), , drop = FALSE]
  ids <- sort(intersect(unique(dat_surv$id), unique(dat_long$id)))
  dat_surv <- dat_surv[match(ids, dat_surv$id), , drop = FALSE]
  dat_long <- dat_long[dat_long$id %in% ids, , drop = FALSE]
  
  X_all <- model.matrix(long_formula, data = dat_long)
  Z_all <- model.matrix(surv_formula, data = dat_surv)
  
  Xs <- s_matrix(X_all, protect = c("obstime"))
  Zs <- s_matrix(Z_all, protect = NULL)
  
  idx_by_id <- split(seq_len(nrow(dat_long)), dat_long$id)
  Xlist <- vector("list", length(ids))
  ylist <- vector("list", length(ids))
  w_index <- vector("list", length(ids))
  
  for (i in seq_along(ids)) {
    idx <- idx_by_id[[as.character(ids[i])]]
    Xlist[[i]] <- Xs$M[idx, , drop = FALSE]
    ylist[[i]] <- dat_long$yij[idx]
    w_index[[i]] <- idx
  }
  
  list(
    ids = ids,
    dat_surv = dat_surv,
    dat_long = dat_long,
    X_all = Xs$M,
    Z = Zs$M,
    Xlist = Xlist,
    ylist = ylist,
    idx_by_id = w_index,
    time = pmax(dat_surv$si, 1e-12),
    event = dat_surv$event,
    X_scale = Xs,
    Z_scale = Zs
  )
}