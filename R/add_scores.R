add_scores <- function(long, surv, grid_cols_X, grid_cols_W,
                       gp, rfpca_nbasis = 15, rfpca_npc_max = NULL) {
  
  Xmat_long <- as.matrix(long[, grid_cols_X, drop = FALSE])
  Wmat_surv <- as.matrix(surv[, grid_cols_W, drop = FALSE])
  
  FPCAp_X <- process_rfpca(long$ID, Xmat_long, gp = gp, 
                           nbasis = rfpca_nbasis, npc_max = rfpca_npc_max)
  FPCAp_W <- process_rfpca(surv$ID, Wmat_surv, gp = gp, 
                           nbasis = rfpca_nbasis, npc_max = rfpca_npc_max)
  
  long_p <- long
  surv_p <- surv
  long_p$FPCscore.X <- I(FPCAp_X$FPCscore)
  surv_p$FPCscore.W <- I(FPCAp_W$FPCscore)
  
  list(
    robust    = list(long = long_p, surv = surv_p, FPCA.X = FPCAp_X, FPCA.W = FPCAp_W)
  )
}
