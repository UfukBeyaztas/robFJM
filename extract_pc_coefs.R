extract_pc_coefs <- function(coef_vec, prefix = "FPCscore.X") {
  nm <- names(coef_vec)
  if (is.null(nm)) return(numeric(0))
  
  idx <- grep(paste0("^", prefix), nm)
  if (length(idx) == 0) return(numeric(0))
  
  subnm <- nm[idx]
  
  dig <- suppressWarnings(as.integer(gsub(".*?(\\d+)$", "\\1", subnm)))
  ord <- if (all(is.finite(dig))) order(dig) else seq_along(idx)
  as.numeric(coef_vec[idx][ord])
}