# Extract group-level fixed effect estimates using LASSO
# 
# Parameters:
#   X               A data frame or matrix with all of the non-fixed controls
#   Y               A vector of the outcome variable
#   W               A vector of the treatment variable
#   FE              A data frame, matrix, or vector containing identifiers for the groups 
#                     for which fixed effects will be extracted
#   other.FE        A data frame, matrix, or vector containing identifiers for fixed controls 
#                     for which fixed effects will not be extracted. Optional
#   degree          The degree of polynomial for the controls
#   orthogonalize   If TRUE, the control variables will be orthogonalized before computing
#                     the polynomials

fe.extract <- function(X, Y, W, FE, other.FE = NULL, degree = 3, orthogonalize = FALSE) {
	
	# load required libraries
	require(Matrix)
	require(glmnet)
  require(dplyr)

  if (missing(X) | missing(Y) | missing(W) | missing(FE)) {
    stop("Please provide values for X, Y, W, and FE.")
  }
  if (degree < 1) {
    stop("degree must be a positive integer.")
  }

  X <- as.data.frame(X)
  Y <- data.matrix(Y)
  W <- data.matrix(W)
  
	# what level of polynomial? why stop at 3? can we make it flexible?
  xmat <- data.matrix(data.frame(do.call(polym, c(X, degree = degree, raw = !orthogonalize))))

	# fixed effects to include but not keep
  ofemat <- NULL
  if (!is.null(other.FE)) {
    if (is.vector(other.FE)) {
      ofemat <- t(as(factor(other.FE), "sparseMatrix"))
    } else {
      for (f in 1:dim(other.FE)[2]) {
        ofemat <- cbind(ofemat, t(as(factor(other.FE[,f]), "sparseMatrix")))
      }
    }
  }

	# create dummy variables for the fixed effects we want to keep
	num.fe <- NULL
	fe.names <- NULL
	femat <- NULL
  if (is.vector(FE)) {
    fe.names <- "fe"
    femat <- t(as(factor(FE), "sparseMatrix"))
    num.fe <- length(unique(FE))
    colnames(femat) <- unique(FE)
  } else {
    fe.names <- colnames(FE)
    for (f in 1:dim(FE)[2]) {
      tempmat <- t(as(factor(FE[,f]), "sparseMatrix"))
      colnames(tempmat) <- paste0(fe.names[f],unique(FE[,f]))
      num.fe <- c(num.fe, length(unique(FE[,f])))
      femat <- cbind(femat, tempmat)
    }
  }

  # combine above into the data
  x_lasso <- Matrix(data.matrix(cbind(W,xmat,ofemat,femat)), sparse = TRUE)
  y_lasso <- Matrix(Y, sparse = TRUE)
  n.not.id = ncol(x_lasso) - sum(num.fe)
  
  # run the LASSO (include OLS version too and make it a switch)
  out.lasso <- cv.glmnet(x=x_lasso,y=y_lasso,alpha=1,penalty.factor=c(rep(1,n.not.id),rep(0,sum(num.fe))))
  
  # pull out the fixed effect estimates we care about
  if (is.vector(FE)) {
    fe.coeffs <- coef(out.lasso, s = "lambda.min")[(n.not.id+2):(sum(num.fe)+n.not.id+1),]
    old.df <- data.frame(fe = as.character(FE))
    new.df <- data.frame(fe = as.character(names(fe.coeffs)),
                         fe.est = fe.coeffs)
    fe.out <- left_join(x = old.df, y = new.df, by = "fe")["fe.est"]
  } else {
    fe.out <- data.frame(matrix(0,dim(Y)[1],length(num.fe)))
    out.names <- paste0(fe.names,".est")
    names(fe.out) <- out.names
    all.coeffs <- coef(out.lasso, s = "lambda.min")[(n.not.id+2):(sum(num.fe)+n.not.id+1),]
    fe.coeffs <- all.coeffs[1:num.fe[1]]
    old.df <- data.frame(fe = as.character(FE[,1]))
    new.df <- data.frame(fe = as.character(sub(fe.names[1],"",names(fe.coeffs))),
                         fe.est = fe.coeffs)
    names(new.df) <- c("fe",out.names[1])
    fe.out[out.names[1]] <- left_join(x = old.df, y = new.df, by = "fe")[[out.names[1]]]
    for (f in 2:length(num.fe)) {
      fe.coeffs <- all.coeffs[(1+sum(num.fe[1:(f-1)])):sum(num.fe[1:f])]
      old.df <- data.frame(fe = as.character(FE[,f]))
      new.df <- data.frame(fe = as.character(sub(fe.names[f],"",names(fe.coeffs))),
                           fe.est = fe.coeffs)
      names(new.df) <- c("fe",out.names[f])
      fe.out[out.names[f]] <- left_join(x = old.df, y = new.df, by = "fe")[[out.names[f]]]
    }
  }

	fe.out
	
}
