"Hello World This is gene mapping project"
"Pleas"
dat <- matrix( 0, 8, 4 )
colnames( dat ) <- c("CCStatus", "XRAY", "BINYOB", "weights")
dat[5:8,1] <- 1
dat[c(3,4,7,8),2] <- 1
dat[c(2,4,6,8),3] <- 1