#Additive model(Default in plink)
A = read.table("plink.mds",header = TRUE) # Read the mds file from plink output
plot(A$C1,A$C2)
CasesPC1<-c(A$C1[1:846])
CasesPC2<-c(A$C2[1:846])
plot(CasesPC1, CasesPC2)
ControlsPC1<-c(A$C1[847:1205])
ControlsPC2<-c(A$C2[847:1205])
points(ControlsPC1,ControlsPC2,col = "red")

B = read.table("plink.assoc.logistic",header = TRUE)
plot(B$BP,-log10(B$P))
library(gap) # Create the Manhattan Plot(gap version)
ord <- with(d,order(B$CHR,B$BP))
test <- data.frame(chr = B$CHR, pos = B$BP, p = B$P)
mhtplot(test)
library(qqman)
op<-(ylim=c(0,10))
manhattan(P,ylim=c(0,10)) # Create the Manhattan Plot(qqman version)
qq(P$P)
OrigGenData <-B
P<-OrigGenData[OrigGenData$TEST == "ADD" & OrigGenData$P,]
# plink --bfile CAMPdata2 --logistic --sex --covar pheno.txt --covar-name GENDER --adjust
# --covar plink.mds --covar-name C1 
nrows = nrow(OrigGenData)


GenDat<-data.frame(chr = CHR, pos = BP, p = P)