###### NOTE ####
## early work in progress
## no guarantees of usefulness
######



# this script correlates the expression values obtained from qPCR and microarray on the same samples
# where your expression values are in a dataframe by sample, including values for log2 microarray and log2 qPCR     
# and currently only functions on three hypothetical genes (to be fixed), sg1, pom210 and setd7

# build a dummy dataset
temp <- rnorm(n = 10, mean = 2, sd = 0.5)
temp2 <- matrix(c(
  temp, temp+rnorm(n=10,mean=0.2,sd=0.5), temp,
  temp+5, sqrt(temp), temp^2), 
          nrow = 10, ncol = 6)
dim(temp2)
colnames(temp2) <- c("C001R002","C001R003","C001R004", "sg1","setd7","pom210")
data <- as.data.frame(temp2)
data

plot(data[,2], data[,5])

### Important Assumptions ###
# currently assumes data is contained in a single dataframe, with the
# order of the probes followed by the qpcr values in the same order 
# i.e. three probes, then three qPCR values in the same gene of interest order 

# here are the names of your probes and qpcr values 
probes <- colnames(data)[1:3]
genes <- colnames(data)[4:6]

cor(x = data[,1], y = data[,5], 
    method = "pearson")


# identify correlation via cor()
num.comps <- length(colnames(data))/2
cor.result <- rep(0, times = num.comps)
for(i in 1:3) {
  cor.result[i] <- cor(x = data[,i], y = data[,i+3])
  names(cor.result)[i] = genes[i]
}

cor.result

#### in order to easily obtain slope, could do this with linear model instead of
# cor, given that the two variables have similar variance

num.comps <- length(colnames(data))/2
lm.result <- list()
for(i in 1:3) {
  lm.result[[i]] <- lm(data[,i] ~ data[,i+3])
  names(lm.result)[i] = genes[i]
}

r = NULL
slope = NULL
for(i in 1:3) {
  r[i] <- sqrt((summary(lm.result[[i]])$r.squared))
  slope[i] <- summary(lm.result[[i]])$coefficients[2]
  names(r)[i] = paste(genes[i],".r", sep = "")
  names(slope)[i] = paste(genes[i],".slope", sep = "")
}

# results:
r
slope

# to confirm the slope, plot the data
par(mfrow = c(2,2))
for(i in 1:3) {
  plot(data[,i+3], data[,i], xlim = c(1,10), ylim = c(1,10))
}

# interleave the data
library(ggplot2)
summary.plot.data <- ggplot2:::interleave(r,slope)

# construct barplot of r value + slope
#par(mfrow = c(1,1))
#generate summary figure for correlations
barplot(summary.plot.data, 
        #beside = TRUE, 
        col=c("black","grey"),
        las = 3, 
        #yaxt = "n"
        )
#axis(2, at = seq(from=0, to=1.4, by= 0.25), las =1)
#legend(20,y=1.4, legend = c("Rsq", "Slope"),
       fill = c("black","grey"), cex = 0.8)