library(sva)

expr <- read.table("finalMerged.txt",row.names="external_gene_id",header=T)
expr2 <- expr[complete.cases(expr),]

### pca before ComBat
pca <- prcomp(t(expr2), retx=T, center=T, scale=T)
pc <- pca$x

tissue <- as.character(c(rep("00",36),sapply(names(expr2)[37:283],function(x){substr(x,14,15)})))
subBatch <- as.character(c(rep("cellLine",36),sapply(names(expr2)[37:283],function(x){substr(x,22,25)})))

plot(pc[,1], pc[,2], col=as.factor(tissue))
plot(pc[,1], pc[,2], col=as.factor(subBatch))

### ComBat
batch <- c(rep("cellLine",36),rep("tcga",247))
expr_combat <- ComBat(dat=expr2, batch, mod=NULL,numCovs=NULL, par.prior=TRUE,prior.plots=FALSE)

### pca after ComBat
pca <- prcomp(t(expr_combat), retx=T, center=T, scale=T)
pc <- pca$x

tissue <- as.character(c(rep("00",36),sapply(names(expr_combat)[37:283],function(x){substr(x,14,15)})))
subBatch <- as.character(c(rep("cellLine",36),sapply(names(expr_combat)[37:283],function(x){substr(x,22,25)})))

plot(pc[,1], pc[,2], col=as.factor(tissue))
plot(pc[,1], pc[,2], col=as.factor(subBatch))

# ### subset cell line data to the vpa treatment and ctr using the expr_combat
# cellline <- expr_combat[,1:36]
# vpa.cellline <- sort(names(cellline)[c(grep("vpa",names(cellline)),grep("ctr",names(cellline)))])[-c(3,4,11,12,19)]
# vpa <- cellline[,vpa.cellline]
# 
# vpa.matrix <- as.matrix(vpa)
# group <- rep(c(1,0),length=16)
# pair <- as.factor(rep(1:8,each=2,length=16))
# 
# #ctr.idx <- seq(1,16,by=2)
# #vpa.idx <- seq(2,16,by=2)
# #t.test(vpa.matrix[1,ctr.idx],vpa.matrix[1,vpa.idx],paired=T)
# 
# ### select significant genes using expr_combat
# result <- matrix(nrow=nrow(vpa.matrix),ncol=4)
# colnames(result) <- c("est_beta0","est_beta1","tstat", "pvalue")
# rownames(result) <- row.names(vpa)
# for (i in 1:nrow(vpa.matrix)){
#   #if(i%%100==0){print(i)}
#   lm1 <- lm(vpa.matrix[i,] ~ group + pair)
#   result[i,] <- c(summary(lm1)$coef[1,1],summary(lm1)$coef[2,c(1,3,4)])
# }
# 
# topGenes <- result[order(result[,4]),]
# fdr <- p.adjust(topGenes[,4],method="fdr")
# topGenes <- cbind(topGenes,fdr)

### subset cell line data to the vpa treatment and ctr using the expr2
cellline <- expr2[,1:36]
vpa.cellline <- sort(names(cellline)[c(grep("vpa",names(cellline)),grep("ctr",names(cellline)))])[-c(3,4,11,12,19)]
vpa <- cellline[,vpa.cellline]

vpa.matrix <- as.matrix(vpa)
group <- rep(c(1,0),length=16)
pair <- as.factor(rep(1:8,each=2,length=16))

#ctr.idx <- seq(1,16,by=2)
#vpa.idx <- seq(2,16,by=2)
#t.test(vpa.matrix[1,ctr.idx],vpa.matrix[1,vpa.idx],paired=T)

# ### select significant genes using expr2
# result <- matrix(nrow=nrow(vpa.matrix),ncol=4)
# colnames(result) <- c("est_beta0","est_beta1","tstat", "pvalue")
# rownames(result) <- row.names(vpa)
# for (i in 1:nrow(vpa.matrix)){
#   #if(i%%100==0){print(i)}
#   lm1 <- lm(vpa.matrix[i,] ~ group + pair)
#   result[i,] <- c(summary(lm1)$coef[1,1],summary(lm1)$coef[2,c(1,3,4)])
# }
# 
# topGenes <- result[order(result[,4]),]
# fdr <- p.adjust(topGenes[,4],method="fdr")
# topGenes <- cbind(topGenes,fdr)

# alternative limma
library(limma)

design <-  model.matrix(~group+pair)
fit <- lmFit(vpa.matrix,design)
fit <- eBayes(fit)

topGenes_vpa <-topTable(fit,coef=2,number=nrow(vpa.matrix))
### the gene selected from the combat adjusted data is the same as the combat unadjusted data (est_beta is different, but the t-statistic and p-values are the same.)

### ASSIGN
library(ASSIGN)
nTop <- 200
geneList <- rownames(topGenes)[1:nTop]
#testData_sub <- expr_combat[geneList,37:283]
testData_sub <- expr2[geneList,37:283]
B_vector <- topGenes[1:nTop,1]
S_matrix <- topGenes[1:nTop,2]
Pi_matrix <- rep(0.95,nTop)

##limma
geneList <- topGenes_vpa$ID[1:nTop]
#testData_sub <- expr_bombat[geneList,37:283]
testData_sub <- expr2[geneList,37:283]
B_vector <- fit$coef[geneList,1]
S_matrix <- fit$coef[geneList,2]
Pi_matrix <- rep(0.95,nTop)


#test1: adaptive_B=F, adaptive_S=F, mixture_beta=F
mcmc.chain <- assign.mcmc(Y = testData_sub, Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=F, adaptive_S=F, mixture_beta=F, p_beta = 0.5)
mcmc.pos.mean1 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=F, adaptive_S=F,mixture_beta=F)

#test2: adaptive_B=T, adaptive_S=F, mixture_beta=F
mcmc.chain <- assign.mcmc(Y = testData_sub, Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=T, adaptive_S=F, mixture_beta=F, p_beta = 0.5)
mcmc.pos.mean2 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=T, adaptive_S=F,mixture_beta=F)

#test3: adaptive_B=T, adaptive_S=T, mixture_beta=F
mcmc.chain <- assign.mcmc(Y = testData_sub, Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=T, adaptive_S=T, mixture_beta=F, p_beta = 0.5)
mcmc.pos.mean3 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=T, adaptive_S=T,mixture_beta=F)

#test4: adaptive_B=T, adaptive_S=T, mixture_beta=T
mcmc.chain <- assign.mcmc(Y = testData_sub, Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=T, adaptive_S=T, mixture_beta=T, p_beta = 0.5)
mcmc.pos.mean4 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=T, adaptive_S=T,mixture_beta=T)

# plots and tables
label <- as.factor(c(rep("normal",32),rep("tumor",215)))
pdf("ComBat_test.pdf")
par(mfrow=c(2,2))
boxplot(mcmc.pos.mean1$beta ~ label,ylab="vpa signature",main="test1")
boxplot(mcmc.pos.mean2$beta ~ label,ylab="vpa signature",main="test2")
boxplot(mcmc.pos.mean3$beta ~ label,ylab="vpa signature",main="test3")
boxplot(mcmc.pos.mean4$kappa ~ label,ylab="vpa signature",main="test4")
dev.off()

pa <- matrix(mcmc.pos.mean4$kappa,nrow=ncol(testData_sub),1);rownames(pa) <- names(testData_sub)
pa1 <- pa[order(rownames(pa)),1,drop=F]
pa2 <- pa1[c(grep("11A",rownames(pa1))-1, grep("11B",rownames(pa1))-1,grep("11A",rownames(pa1)), grep("11B",rownames(pa1))),1,drop=F]
dim(pa2) <- c(32,2)
rownames(pa2) <- sapply(rownames(pa1)[1:32],function(x){substr(x,1,12)})
colnames(pa2) <- c("tumor", "normal")

write.csv(pa2,file="Combat_test4.csv")

