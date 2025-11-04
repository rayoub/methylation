function (x, batch)
{
	batch <- as.factor(batch)
	contrasts(batch) <- contr.sum(levels(batch))
	batch <- model.matrix(~batch)[, -1, drop = FALSE]
    X.batch <- batch
	design <- matrix(1, ncol(x), 1)
    x <- as.matrix(x)
    fit <- lmFit(x, cbind(design, X.batch), ...)
    beta <- fit$coefficients[, -(1:ncol(design)), drop = FALSE]
    beta[is.na(beta)] <- 0
    x - beta %*% t(X.batch)
}