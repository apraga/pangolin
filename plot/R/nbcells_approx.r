x <- seq(1, 90)

dlat <- 1
coef <- pi/180.

y_approx <- c()
y_trunc <- c()
y_true <- c()
for (i in seq(1, 90)) {
  y_approx[i] <- 3*(2*i-1)
  y_trunc[i] <- min(3*(2*i-1), 3*(2*60-1))
  y_true[i] <- 3*(2*sin(0.5*dlat*coef)*sin((i*dlat)*coef))/(1-cos(dlat*coef))
}

postscript("nbcells_approx.ps",width=9,height=9,pointsize=15)

# Revert sorting
ylim <- c(min(y_approx), max(y_approx))
plot(x, y_approx, type='l', ylim=ylim, xaxt='n', xlab='', ylab='', col='red', 
    lwd=2)
par(new=TRUE)
plot(x, y_trunc, type='l', ylim=ylim,xaxt='n', xlab='', ylab='', col='blue', 
    lwd=2, lty=2)
par(new=TRUE)
plot(x, y_true, type='l', ylim=ylim, xaxt='n', xlab='latitude', ylab='nb cells', 
    lwd=2, col='black', lty=3)
xlab <- seq(90, 0, -10)
axis(1, at=rev(xlab), label=xlab)

legend_pos = "topleft"
legend(legend_pos, c('Pangolin', 'Truncated', 'Exact'), col=c('red', 'blue', 'black'), 
    lty=c(1,2,3), lwd=2)


title("Number of cells approximation (1x0.67°)")
dev.off()
