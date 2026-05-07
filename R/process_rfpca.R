process_rfpca <- function(subjs, funcs_mat, gp, nbasis = 15, npc_max = NULL) {
  ids <- as.numeric(as.factor(subjs))
  base_idx <- !duplicated(ids)
  baseline <- as.matrix(funcs_mat[base_idx, , drop = FALSE])
  
  rp <- getRPCA(baseline, nbasis = nbasis, gp = gp)
  
  if(is.null(npc_max)){
    K <- rp$ncomp
  }else{
    K <- min(rp$ncomp, npc_max)
  }
  
  if (!is.finite(K) || K < 1) K <- 1
  
  scores0 <- rp$PCAscore[, 1:K, drop = FALSE]
  colnames(scores0) <- paste0("PC", 1:K)
  
  scores <- scores0[ids, , drop = FALSE]
  colnames(scores) <- paste0("PC", 1:K)
  
  basis_eval <- eval.fd(gp, rp$PCAcoef)[, 1:K, drop = FALSE]
  colnames(basis_eval) <- paste0("PC", 1:K)
  
  list(FPCscore = scores, kz = K, basis = basis_eval, rp = rp)
}
