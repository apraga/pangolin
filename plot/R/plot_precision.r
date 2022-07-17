# Plot precision for nb_lat2=90 for a varying CFL

dir <- "/wkdir/pae2/praga/precision"


names <- c('0_9', '0_75', '0_7', '0_5', '0_25', '0_1', 
           '432', '518', '555', '777', '1555', '3887')


names <- matrix(names, ncol=2)

# Subset
nstart <- 1
nend <- dim(names)[1]

mse <- vector()
mse_nolimit <- vector()
mse_nolimit2 <- vector()
# Mean squared error
MSE <- function(num) mean(num^2); 


x <- seq(nstart, nend)

for (i in x) {
    name <- paste(dir, "/ratio_anal_CFL_", names[i,1], "_", names[i,2], ".dat", sep="")
    ref <- read.table(name)

    name <- paste(dir, "/ratio_seq_CFL_", names[i,1], "_", names[i,2], ".dat", sep="")
    cur <- read.table(name)
    diff <- ref[1] - cur[1]
    mse[i] <- MSE(diff)

    name <- paste(dir, "/ratio_seq_nolimit_CFL_", names[i,1], "_", names[i,2], ".dat", sep="")
    cur_nolimit <- read.table(name)
    diff <- ref[1] - cur_nolimit[1]
    mse_nolimit[i] <- MSE(diff)
}

CFL <- c(0.9, 0.75, 0.7, 0.5, 0.25, 0.1)
#postscript("precision_cfl_both.ps",width=10,height=10,pointsize=15)
lim <-c (min(min(mse), min(mse_nolimit)), max(max(mse), max(mse_nolimit))) 
plot(CFL, mse, type="o", ylab="Mean squared error", ylim=lim)
par(new=TRUE)
plot(CFL, mse_nolimit, type="o", lty=2, ylab="Mean squared error", ylim=lim)
legend("topright", c("with slope limitation", "without"), lty=c(1,2))
title("Precision for a varying timestep")
#dev.off()
