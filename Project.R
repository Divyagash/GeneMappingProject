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

# Adjusted for the first PC
AdjGenData = read.table("plink.assoc.logistic",header = TRUE)
op<-(ylim=c(0,10))
manhattan(P,ylim=c(0,10))
qq(P$P)
OrigGenData <-B
NewAdj<-AdjGenData[AdjGenData$TEST == "C1" & AdjGenData$P,]
NewAdj2<-NewAdj[complete.cases(NewAdj),]
op<-(ylim=c(0,10))
manhattan(NewAdj2,ylim=c(0,10),main="Manhattan plot after adjusting for the first PC")
qq(NewAdj2$P)

# Outlier analysis
Outlier = read.table("plink.nearest",header = TRUE)
head(Outlier)
hist(Outlier$Z)
min(Outlier$Z)
NN1<-Outlier[Outlier$NN == 1 & Outlier$Z,]
hist(NN1$Z, xlab = "Z-score", main = "Histogram of First Nearest Neighbor Z-scores")
abline(v=-3,col ="red")
Outlier_1stNN<-NN1$IID[NN1$Z < -3]
Outlier_1stNNFD<-NN1$FID[NN1$Z < -3]

NN2<-Outlier[Outlier$NN == 2 & Outlier$Z,]
Outlier_2ndNN<-NN2$IID[NN2$Z < -3]
Outlier_2ndNNFD<-NN2$FID[NN2$Z < -3]

hist(NN2$Z,xlab = "Z-score", main = "Histogram of Second Nearest Neighbor Z-scores")
abline(v=-3,col ="red")
NN3<-Outlier[Outlier$NN == 3 & Outlier$Z,]
Outlier_3rdNN<-NN3$IID[NN3$Z < -3]
Outlier_3rdNNFD<-NN3$FID[NN3$Z < -3]

hist(NN3$Z,xlab = "Z-score", main = "Histogram of Third Nearest Neighbor Z-scores")
abline(v=-3,col ="red")
NN4<-Outlier[Outlier$NN == 4 & Outlier$Z,]
Outlier_4thNN<-NN4$IID[NN4$Z < -3]
Outlier_4thNNFD<-NN4$FID[NN4$Z < -3]

hist(NN4$Z,xlab = "Z-score", main = "Histogram of Fourth Nearest Neighbor Z-scores")
abline(v=-3,col ="red")

NN5<-Outlier[Outlier$NN == 5 & Outlier$Z,]
Outlier_5thNN<-NN5$IID[NN5$Z < -3]
Outlier_5thNNFD<-NN5$FID[NN5$Z < -3]

hist(NN5$Z,xlab = "Z-score", main = "Histogram of Fifth Nearest Neighbor Z-scores")
abline(v=-3,col ="red")


Subj1<-Outlier[Outlier$IID == "SUBJ-0015" & Outlier$Z,]
hist(Subj1$Z)



#Related Analysis
Related = read.table("plink.genome",header = TRUE)
Pi_hat <- Related$PI_HAT
abline(v=0.05)
max(Pi_hat)
hist(Pi_hat,main = "Histogram of Proportion of IBD", xlab = "IBD proportion")
Rel<-Related$IID1[Pi_hat > 0.1]
Rel2<-Related$IID2[Pi_hat > 0.05]
hist(Related$PPC)
related <-c()

for (i in 1:length(Rel)){
  if ((Rel[i]  %in% related) == FALSE){
    related<-c(related, Rel[i])
    
  }
  
}
unique(unlist(related))
A<-table(related)
A[A>100]

relatedDF<-cbind(Rel,Rel2)

CoupRelated<-data.frame()
for (i in 1:length(Rel)){
  coup<-data.frame(Rel = Rel2[i], Rel2 = Rel[i])
  Bool<-c(TRUE, TRUE)
  if (nrow(merge(relatedDF,coup))>0){
    CoupRelated<-rbind(CoupRelated,coup)
    
    
    
  } 
}
CoupRelated
#SUBJ-0126 SUBJ-0209 SUBJ-0420 SUBJ-0514 SUBJ-0617 SUBJ-1082

#Plot depicting the relationship between the subject
Z0<-Related$Z0
Z1<-Related$Z1
plot(Z0,Z1,xlim=c(0,1),ylim=c(0,1), main = "Relatedness among subjects")
abline(1,-1)

