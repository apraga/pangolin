dir <- "/wkdir/pae2/praga/hourdin"
#dir <- "output"
ref <- list()
ref <- c(ref, list(read.table(paste(dir, "/ratio_anal_220_25.dat", sep=""))))
ref <- c(ref, list(read.table(paste(dir, "/ratio_anal_263_30.dat", sep=""))))
ref <- c(ref, list(read.table(paste(dir, "/ratio_anal_392_45.dat", sep=""))))
ref <- c(ref, list(read.table(paste(dir, "/ratio_anal_777_90.dat", sep=""))))
ref <- c(ref, list(read.table(paste(dir, "/ratio_anal_1163_135.dat", sep=""))))
#ref <- c(ref, list(read.table(paste(dir, "/ratio_anal_1549_180.dat", sep=""))))
#ref <- c(ref, list(read.table(paste(dir, "/ratio_anal_1935_225.dat", sep=""))))

q <- list()
q <- c(q, list(read.table(paste(dir, "/tmp25/ratio_seq_220.dat", sep=""))))
q <- c(q, list(read.table(paste(dir, "/tmp30/ratio_seq_263.dat", sep=""))))
q <- c(q, list(read.table(paste(dir, "/tmp45/ratio_seq_392.dat", sep=""))))
q <- c(q, list(read.table(paste(dir, "/tmp90/ratio_seq_777.dat", sep=""))))
q <- c(q, list(read.table(paste(dir, "/tmp135/ratio_seq_1163.dat", sep=""))))
#q <- c(q, list(read.table(paste(dir, "/tmp180/ratio_seq_1549.dat", sep=""))))
#q <- c(q, list(read.table(paste(dir, "/tmp225/ratio_seq_1935.dat", sep=""))))

diff <- list()
# Append to the end
for (i in seq(1,length(ref))) {
diff <- c(diff, list(ref[[i]][1] - q[[i]][1]))
}

qmse <- vector()
# Mean squared error
MSE <- function(num) mean(num^2); 

for (i in seq(1,length(ref))) {
## 2 brackets
    #qmean[i] <- abs(colMeans(diff[[i]]))
    qmse[i] <- MSE(diff[[i]])
}
#
## Custom x axis
xlabel <- c('3.6x2.45','3x2.03','2x1.34','1x0.67','0.67x0.45')#, '0.5x0.33')#, '0.4x0.27')
x <- c(3.6*2.45,3*2.03,2*1.34,1*0.67,0.67*0.45)#,0.5*0.33)#,0.4*0.27)
x <- rev(x)
#x <- seq(1, length(ref))
#xlabel <- x

png("mse.png",width=10,height=10,units="cm",res=600,pointsize=8)
plot(x,qmse,type='o',xaxt='n', ylab='')
axis(side=1,at=x,labels=xlabel)
title(ylab="MSE")
dev.off()
