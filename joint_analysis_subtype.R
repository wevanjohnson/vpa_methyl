## integrate methylation and expression results by subtypes
met <- read.csv("methylation_tcga_vpa_all_subtype.csv",as.is=T)
sampleID <- as.character(c(met[1:32,1],met[1:213,4]))
vpa.m_met <- as.numeric(c(met[1:32,2],met[1:213,5]))
subtype <- as.character(c(met[1:32,3],met[1:213,6]))

met_vpa <- data.frame(sampleID, vpa.m_met, subtype)

expr <- read.table("subTypeSigs.txt",header=T)
expr_vpa <- expr[,c(1,2,3,6)] 
names(expr_vpa) <- c("sampleID", "subtype","vpa.m_expr", "vpa.s_expr")
expr_vpa[,2] <- gsub("Control", "Normal",expr_vpa[,2])
expr_met <- merge(expr_vpa, met_vpa, by=c("sampleID","subtype"))

expr_met <- expr_met[order(expr_met[,2]),]
expr_met$idx <- 1:246

expr_met <- expr_met[c(211:241,1:210,242:246),]
expr_met$idx <- 1:246
expr_met$vpa_expr_met_s <- (expr_met$vpa.s_expr+expr_met$vpa.m_met)/2
expr_met$vpa_expr_met_m <- (expr_met$vpa.m_expr+expr_met$vpa.m_met)/2
expr_met <- expr_met[,-6]
write.csv(expr_met, file="VPA_expr_met_joint.csv")

## methylation signature strength vs subtypes
met <- read.csv("methylation_tcga_vpa_all_subtype.csv",as.is=T)

pdf("methylation_subtype.pdf")


boxplot(met[,c(2,5,8,11,14,17,20)], ylab="ASSIGN Signature", las=1, at=c(1,2,4,5,6,7,8), 
        cex.lab=0.8, col="brown",frame=F, xaxt = "n", main="Methylation Profile")
legend("topleft",
       legend=c("VPA"),
       col=c("brown") ,
       pch=15)
axis(1.,at=c(1,2,4,5,6,7,8),labels=c("Normal","Tumor","Basal","Her2","Lum-A","Lum-B","Normal-Like"),tick=FALSE,cex.axis=0.7)
dev.off()


## (expression + methylation)/2 signature strength vs subtypes

expr_met <- read.csv("VPA_expr_met_joint.csv",as.is=T)[,c(2,3,7,8)]
vpa_s <- expr_met[,3]
#vpa_m <- expr_met[,4]
pdf("Expr_methyl_subtype_m.pdf")
boxplot(list(vpa_s[1:31],vpa_s[32:246],vpa_s[32:70],vpa_s[71:84],vpa_s[85:192],vpa_s[193:241],vpa_s[242:246]), ylab="ASSIGN Signature", las=1, at=c(1,2,4,5,6,7,8), 
        cex.lab=0.8, col="brown",frame=F, xaxt = "n", main="Expression+Methylation Profile")
legend("topleft",
       legend=c("VPA"),
       col=c("brown") ,
       pch=15)
axis(1.,at=c(1,2,4,5,6,7,8),labels=c("Normal","Tumor","Basal","Her2","Lum-A","Lum-B","Normal-Like"),tick=FALSE,cex.axis=0.7)
dev.off()




#####################
expr_met <- read.csv("VPA_expr_met_joint.csv",as.is=T)

pdf("boxplot_subtype_joint_s.pdf",height=7, width=10)
boxplot(expr_met$vpa_expr_met_s ~ factor(expr_met$subtype, levels=c("Normal","Basal","Her2","LumA","LumB","Normal-Like")), ylab="ASSIGN Signature", las=1,  
        names=c("Normal","Basal","Her2","Lum-A","Lum-B","Normal-Like"),
        cex.lab=0.8, frame=F)
dev.off()
