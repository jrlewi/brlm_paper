

fit <- nTheory
plot(fit$Beta, pch = 19, cex = .3)
plot(fit$betal[1,], pch = 19, cex = .3)
plot(fit$betal[2,], pch = 19, cex = .3)
plot(fit$Beta,fit$betal[2,], pch = 19, cex = .3)
cor(fit$Beta,fit$betal[2,])
plot(fit$betal[2,],fit$betal[1,])

plot(fit$sigma2s[,1], pch = 19, cex = .3)
plot(fit$sigma2s[,2], pch = 19, cex = .3)
plot(fit$bstar, pch = 19, cex = .3)
plot(fit$mu_rho, pch = 19, cex = .3)
plot(fit$psi_rho, pch = 19, cex = .3)
plot(fit$rho, pch = 19, cex = .3)



plot(marginals[rest_ind, i,], marginals[rlm_ind, i,], pch = 19, cex = .1)
abline(0,1, lwd =2 , col = 2)


plot(marginals[rest_ind, i,], marginals[ols_ind, i,], pch = 19, cex = .1)
abline(0,1, lwd =2 , col = 2)

plot(marginals[rest_ind, i,], marginals[t_ind, i,], pch = 19, cex = .1)
abline(0,1, lwd =2 , col = 2)

predictions[rlm_ind, i, ]

plot(predictions[rest_ind, i, ], predictions[rlm_ind, i, ], pch = 19, cex = .1)
abline(0,1, lwd =2 , col = 2)



plot(group_estimates[rest_ind, 1:p, , i],group_estimates[rlm_ind, 1:p, , i], pch = 19, cex = .5)
abline(0,1, lwd =2 , col = 2)
)

plot(group_estimates[rest_ind, p+1, , i],group_estimates[rlm_ind, p+1, , i], pch = 19, cex = .5)
abline(0,1, lwd =2 , col = 2)
)