getRPCA <- function(data, nbasis, gp){
  n <- dim(data)[1]
  p <- dim(data)[2]
  dimnames(data) = list(as.character(1:n), as.character(1:p))
  bs_basis <- create.bspline.basis(rangeval = c(gp[1], gp[p]), nbasis = nbasis)
  inp_mat <- inprod(bs_basis, bs_basis)
  sinp_mat <- sqrtm(inp_mat)
  evalbase = eval.basis(gp, bs_basis)
  pcaobj <- smooth.basisPar(gp, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  
  mean_coef <- l1median(t(pcaobj$coefs), trace = -1)
  sdata <- scale(t(pcaobj$coefs), center = mean_coef, scale = FALSE)
  new.data <- sdata %*% sinp_mat
  dcov <- covPCAproj(new.data)
  d.eigen <- eigen(dcov$cov)
  var_rob <- cumsum(d.eigen$values)/sum(d.eigen$values)
  ncomp <- which(var_rob > 0.95)[1]
  ppur <- PCAproj(new.data, ncomp)
  loads <- ppur$loadings
  PCs <- solve(sinp_mat) %*% loads
  colnames(PCs) = 1:ncomp
  for(i in 1:ncomp)
    colnames(PCs)[i] = paste("PC", i, sep = "")
  PCAcoef <- fd(PCs, bs_basis)
  mean_coef <- fd(as.vector(mean_coef), bs_basis)
  pcaobj2 <- pcaobj
  pcaobj2$coefs <- t(sdata)
  PCAscore <- inprod(pcaobj2, PCAcoef)
  colnames(PCAscore) = 1:ncomp
  for(i in 1:ncomp)
    colnames(PCAscore)[i] = paste("PC", i, sep = "")
  
  return(list(PCAcoef = PCAcoef, PCAscore = PCAscore, meanScore = mean_coef,
              bs_basis = bs_basis, evalbase = evalbase, ncomp = ncomp, gp = gp))
}