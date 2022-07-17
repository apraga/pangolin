#Plot scheme convergence with mean squared error vs resolution
# Needs the difference to be stored already (see create_diff.r)

dir <- "/wkdir/pae2/praga/convergence"


names <- c('025', '030', '045', '090', '135', '180', '225', '270', '360', '450',
           '675',
           '220', '263', '392', '777', '1163', '1549', '1935', '2320', '3092', 
           '3864', '5792')

names <- matrix(names, ncol=2)

# Subset
nstart <- 1
nend <- dim(names)[1]

qmse <- vector()
# Mean squared error
MSE <- function(num) mean(num^2); 

x <- seq(nstart, nend)
for (i in x) {
    #qmean[i] <- abs(colMeans(diff[[i]]))
    cur <- paste(dir, "/diff_seq_", names[i,1], "_", names[i,2], ".dat", sep="")
    diff <- read.table(cur)
    qmse[i] <- MSE(diff)
}
#
## Custom x axis
xlabel <- c('3.6x2.45','3x2.03','2x1.34','1x0.67','0.67x0.45', '0.5x0.33',
    '0.4x0.27','0.33x0.22', '0.25x0.17', '0.2x0.13', '0.13x0.09')
x <- c(3.6*2.45, 3*2.03, 2*1.34, 1*0.67, 0.67*0.45, 0.5*0.33,
    0.4*0.27, 0.33*0.22, 0.25*0.17, 0.2*0.13, 0.13*0.09)
x <- rev(x)

# Avoid to go outside A4 size...
postscript("scheme_error.ps",width=8,height=8,pointsize=18)
# Space for x title
par(mar=c(7,4,3,2)+0.1, lwd=2)
plot(x[nstart: nend], qmse, type='o', xaxt='n', ylab='', xlab='', log='x',
    lwd=3)
axis(side=1,at=x,labels=xlabel,las=2)
title(ylab="Mean squared error")
# Title under
title(xlab="log(equatorial cell area)", line=5)
title("Scheme convergence")
dev.off()
