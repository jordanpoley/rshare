###### NOTE ####
## early work in progress
######

# this script correlates the expression values obtained from qPCR and microarray on the same samples
# where your expression values are in a dataframe by sample, including values for log2 microarray and log2 qPCR     
# and currently only functions on three hypothetical genes (to be fixed), sg1, pom210 and setd7

# models first to obtain R-squared
sg1 <- lm(results$C042R126 ~ results$sg1)
pom210 <- lm(results$C088R114 ~ results$pom210)
setd7 <- lm(results$C263R087 ~ results$setd7)

sg1.adj.r.squared <- summary(sg1)$adj.r.squared
pom210.adj.r.squared <- summary(pom210)$adj.r.squared
setd7.adj.r.squared <- summary(setd7)$adj.r.squared

rsq <- c(sg1.adj.r.squared, pom210.adj.r.squared, setd7.adj.r.squared)
# names(rsq) <- c("sg1", "pom210", "setd7")
rsq

# plot the correlation
ylim.plot = c(-8,4)
xlim.plot = c(-8,4)
xtxt = -5.5
ytxt = 3.7
par(mfrow=c(1,3), mar = c(4,4,1,2))
plot(probes.for.cor.trimmed$C042R126 ~ {median(results$sg1)- results$sg1}, xlab = "ddCt sg1", ylab = "microarray log2(sg1)", las = 1
     , ylim = ylim.plot, xlim = xlim.plot
     )
text(x = xtxt, y = ytxt, labels = paste("R2 =", round(rsq[1], digits = 2)))
plot(probes.for.cor.trimmed$C088R114 ~ {median(results$pom21)- results$pom210}, xlab = "ddCt pom210", ylab = "microarray log2(pom210)", las = 1 
     , ylim = ylim.plot, xlim = xlim.plot
     )
text(x = xtxt, y = ytxt, labels = paste("R2 =", round(rsq[2], digits = 2)))
plot(probes.for.cor.trimmed$C263R087 ~ {median(results$setd7)- results$setd7}, xlab = "ddCt setd7", ylab = "microarray log2(setd7)", las = 1
     , ylim = ylim.plot, xlim = xlim.plot
     )
text(x = xtxt, y = ytxt, labels = paste("R2 =", round(rsq[3], digits = 2)))
