# Plot error evolution vs time
library(ggplot2)
library(scales) # access date format

dir <-"/wkdir/pae2/praga/parallel/error/hourdin_80lat_1tracer"
cur <- paste(dir, "/error_80lat_1tracer.dat", sep="")
err <- read.table(cur)

date <- c()
colnames(err) <- c("date1", "date2", "l2", "loo")
for (i in seq(length(err$date1))){
  date[i] <- paste(err$date1[i], err$date2[i])
}
date <- as.POSIXlt(date)

# Melt data by hand
labels = c( rep("l2", length(err$l2)), rep("loo", length(err$l2)))
err2 <- data.frame(date=c(date, date), error=c(err$l2, err$loo), type=labels)

pdf("error_evol_hourdin_80lat.pdf")
p <- ggplot(err2, aes(x=date, y=error)) + geom_line() 
p <- p + facet_grid(type~.,scales="free")
p <- p + ggtitle("Analytical error for snail test (1.12x0.75°, 3 cores)")
#p <- p + scale_x_continuous(breaks=as.numeric(err2$date))
p <- p + scale_x_datetime(
    #breaks=c(as.Date("2013-01-01 00:00:00"),as.Date("2013-02-01 03:00:00")))
    breaks=c(as.POSIXct("2013-01-01 00:00"),
        as.POSIXct("2013-01-01 12:00"),
        as.POSIXct("2013-01-02 04:00"),
        as.POSIXct("2013-01-02 19:00")),
    labels = date_format("%m-%d %H:%M"))
print(p)
dev.off()
