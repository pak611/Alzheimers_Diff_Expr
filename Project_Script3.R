
# for plink format
#  https://ida.loni.usc.edu/pages/access/geneticData.jsp#63

# set working directory
setwd("~/Dropbox/Ph.D/Bioinformatics/Project/datasets")


library(BiocManager)
# Mac (terminal): xcode-select --install, to update Rcpp
#install.packages("Rcpp")
#BiocManager::install("snpStats") # Yes -- install from source, update all
library(snpStats)


# load in the genotype data
file_prefix <- "ADNI_GO_2_Forward_Bin"
file.bed <- paste(getwd(),"/",file_prefix,".bed",sep="")
file.bim <- paste(getwd(),"/",file_prefix,".bim",sep="")
file.fam <- paste(getwd(),"/",file_prefix,".fam",sep="")

file.bed

# read the data in plink format
ADNI.data <- read.plink(file.bed, file.bim, file.fam, na.strings = ("-9"))
str(ADNI.data)


genotypes <- ADNI.data$genotypes

# 432 subjects with 730525 SNPs
str(genotypes)


# paper genotypes
# rs6546366
# rs2070852 present!
# rs927010
# rs7768046
# rs3755557
# rs60872856 "A genome-wide association study of plasma phosphorylated tau181

# search for existance of specific SNPs
sum(as.numeric(colnames(genotypes)=="rs60872856"))

# returns Calls, Call.rate, Certain.calls, RAF, MAF, P.AA, P.AB, P.BB, z.WHE
geno.stats <- col.summary(genotypes)
colnames(geno.stats) # Call.rate, MAF, z.HWE, P.AB...
maf.keep <- geno.stats$MAF>.01  # minor allele frequencies
na.rm <- is.na(geno.stats$MAF)   # will need to remove these na's
sum(as.integer(maf.keep),na.rm=T)  # keep about 663,594. 


geno.stats

# Setting thresholds
call <- 0.95
minor <- 0.01



# Filter on MAF and call rate
use <- with(geno.stats, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE                # Remove NA's as well

cat(((ncol(genotypes)-sum(use))/sum(as.integer(maf.keep),na.rm=T))*100,"% of SNPs will be removed due to low MAF or call rate.\n") 
#1960222 SNPs will be removed

# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotypes <- genotypes[,use]
geno.stats <- geno.stats[use,]

print(genotypes) # 661,194 SNPs remain


#remove samples with NA phenotypes and merge with genotype data

# read in the phenotype data (PTAU181)
pheno.data <- read.csv("p_tau181_ADNI_GO_2.csv")

pheno.data
colnames(pheno.data)

library(dplyr)
pheno.subcols <- pheno.data %>% select(RID, PLASMAPTAU181,
)

rownames(genotypes)


pheno.subcols[,"RID"]


# calculate the average PTAU181 level for all visits for each subject

j <- 0
i <- 1

# avg.mat2 will contain the averages. rid.mat2 will contain the RID
avg.mat2 <-c()
rid.mat2 <-c()

while (i < 3758) {
  j <- 0
  while (pheno.subcols[i,'RID']==pheno.subcols[i+j,'RID']) {
    j = j+1
    print("equal")
    print(pheno.subcols[i+j,'RID'])
  }
  print("not equal")
  print(pheno.subcols[i+j,'RID'])
  
  average <- mean(pheno.subcols[(i:(i+j-1)),"PLASMAPTAU181"])
  rid <- pheno.subcols[(i+j-1),"RID"]
  i = i + j
  #print(average)
  #print(i)
  #print(j)
  avg.mat2 <- append(avg.mat2, average)
  rid.mat2 <- append(rid.mat2, rid)
}

hist(avg.mat2)

boxplot(avg.mat2,range=0,ylab="raw probe intensity", main="Raw")

avg.mat2


#install.packages("ggstatsplot")
library(ggstatsplot)
ggbetweenstats(avg.mat2,
                outlier.tagging = TRUE)


Q <- quantile(avg.mat2, probs=c(.25, .75), na.rm = FALSE)

iqr <- IQR(avg.mat2)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range﻿

no_outliers<- subset(avg.mat2, avg.mat2 > (Q[1] - 1.5*iqr) & avg.mat2 < (Q[2]+1.5*iqr))

boxplot(no_outliers)

# log-log transform for normal distribution
#install.packages("moments")
library(moments)
skewness(avg.mat2, na.rm = TRUE)

avg.mat3 <-log10(avg.mat2)

boxplot(avg.mat3)

hist(avg.mat3)


# there are 1190 total subjects
str(rid.mat2)

#************************************************************************************
#*******************[ FROM HERE WE GO TO LINKAGE DISEQUILIBRIUM PRUNING]**************

# comb_mat will merge RID with the cooresponding average
comb_mat<- cbind(rid.mat2,avg.mat2)
comb_list <- list(comb_mat)

# intersect the genotype data with the phenotype data
common_samples <- intersect(rownames(genotypes),comb_mat[,1])

# there are 98 common samples
str(common_samples)

# construct a genotype matrix (genotype2) which only includes the common samples
genotypes2 <- genotypes[common_samples,]

# there are now 98 of the original 432 samples
dim(genotypes2)
dim(genotypes)


common_samples
# now make make the combined matrix comprised of the subsamples which
# are at the intersect of the pheno and geno data
comb_mat.subsamples <- as.data.frame(comb_mat) %>% filter(rid.mat2 %in% common_samples)

str(comb_mat.subsamples)

#5 ld prune with SNPRelate and gds file

ld.thresh <- 0.2    # LD cut-off

# Create gds file, required for SNPRelate functions
#library(BiocManager)
#BiocManager::install("SNPRelate")
library(SNPRelate)
# already ran to create gds file


bed.fn <- "ADNI_GO_2_Forward_Bin.bed"
fam.fn <- "ADNI_GO_2_Forward_Bin.fam"
bim.fn <- "ADNI_GO_2_Forward_Bin.bim"

out.gdsfn <- "ADNIGO6.gds"


gds_file <- snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "ADNIGO6.gds")

gdsfile <- openfn.gds("ADNIGO6.gds", readonly = FALSE)

str(gdsfile$root["snp.position"])

snp_pos <- read.gdsn(index.gdsn(gdsfile, "snp.position"))
snp_id <- read.gdsn(index.gdsn(gdsfile, "snp.id"))
snp_chr <- read.gdsn(index.gdsn(gdsfile, "snp.chromosome"))
map.df <- data.frame(chr=snp_chr,snp_id=snp_id,cm=0,bp=snp_pos) 
map.df <- map.df %>% mutate_at("snp_id", as.character) 
map.df[10000,]


map.df

#LD Prune SNPs for IBD analysis
set.seed(1000)

samp_id <- read.gdsn(index.gdsn(gdsfile, "sample.id"))
snp_id <- read.gdsn(index.gdsn(gdsfile, "snp.id"))


snpSUB <- snpgdsLDpruning(gdsfile, ld.threshold = ld.thresh,
                          #sample.id = rownames(genotypes2), # Only analyze the filtered samples
                          snp.id = colnames(genotypes2)) # Only analyze the filtered SNPs


common.stuff <- intersect(snp_id, colnames(genotypes2))

common.stuff
genotypes

colnames(genotypes)

# names of SNPs to keep from ld pruning
snpset <- unlist(snpSUB, use.names=FALSE)

ld.usesnps <- colnames(genotypes2) %in% snpset 

ld.usesnps

# the number of original SNPs
length(ld.usesnps)

# the number left after LD pruning
sum(ld.usesnps)


# the fraction of SNPs used after LD pruning
sum(ld.usesnps)/length(ld.usesnps)



genotypes2 <- genotypes2[,ld.usesnps]

str(genotypes2)



# now we have 98 subjects with 74,376 SNPs
dim(genotypes2)


#************************************************************************************
#*******************[THE ACTUAL ANALYSIS BEGINS HERE]**********************************



#5. gwas qtl with KQ ratio phenotype and ld/maf filtered snps
library(npdr)
genoNum <- as(genotypes2, "numeric")


# use univariate regression of average PTAU181 regressed onto the genotype using linear model
pheno.fit <- npdr::uniReg(outcome=comb_mat.subsamples$avg.mat2, dataset=genoNum, regression.type="lm")


pheno.fit[1:10,]



genotypes3_mat <- matrix(as.numeric(genotypes2),
                         ncol = ncol(genotypes2))

str(genotypes3_mat)

genotypes3_df <- as.data.frame(genotypes3_mat)

str(genotypes3_df)

str(genotypes3_mat)

# combine the genotype and phenotype data
genotypes3_mat2 <- cbind(genotypes3_mat,comb_mat.subsamples$PTAU181)

# name the phenotype column
colnames(genotypes3_mat2)[ncol(genotypes3_mat2)] <- "PTAU"

# convert to dataframe
genotypes3_df_2 <- as.data.frame(genotypes3_mat2)

# attempt cnCV
cncv.genotype <- consensus_nestedCV(train.ds = genotypes3_df_2, 
                                validation.ds =  NULL, 
                                label = "PTAU",
                                method.model = "regression",
                                is.simulated = FALSE,
                                ncv_folds = c(5, 5),
                                param.tune = FALSE,
                                learning_method = "xgbTree", 
                                importance.algorithm = "ReliefFbestK",
                                relief.k.method = "k_half_sigma",     # surf k
                                num_tree = 10,
                                verbose = F)

#* CNCV THINKS THAT GENOTYPE DATA IS THE SAME ACROSS SAMPLES

cat("\n Train Accuracy [",cncv.genotype$cv.acc,"]\n")
cat("\n Validation Accuracy [",cncv.genotype$Validation,"]\n")
cat("\n Selected Features \n [",cncv.genotype$Features,"]\n")  # features ordered by frequency Reduce(intersect,x))
cat("\n Elapsed Time [",cncv.genotype$Elapsed,"]\n")
# cncv random forest model
cncv.model <- cncv.genotype$Train_model # random forest model trained by cncv
# re-apply it to data
#cncv.predict.df <-data.frame(predict=stats::predict(cncv.model, newdata =), true=)
cncv.predict.df



#************************[ TRY USING LASSO AND RIDGE]************************


x <- genotypes3_mat
y <- comb_mat.subsamples$avg.mat2

comb_mat.subsamples

n<-nrow(genotypes3_mat)
#randomly selects from the entire sets 2/3 of the data
train_rows <- sample(1:n, .66*n)
train_rows
x.train <- x[train_rows, ]
x.test <- x[-train_rows, ]
y.train <- y[train_rows]
y.test <- y[-train_rows]

x.train
y.train

library(glmnet)

# glmnet with alpha=0 for ridge regression

alpha0.fit <-cv.glmnet(x.train, y.train, 
                       type.measure="mse", alpha=0, family="gaussian")


alpha0.fit

alpha0.predicted <- predict(alpha0.fit, s=alpha0.fit$lambda.min, newx=x.test)

alpha0.fit


alpha0.mse <- round(mean((y.test-alpha0.predicted)^2),digits=4)


plot(y.test,alpha0.predicted,
     main="concensus cross validation ridge-regression",
     ylab="actual holdup",xlab="predicted holdup", lines(x,x, col='red'))

# now lets apply elastic net regression (LASSO)

alpha1.fit <-cv.glmnet(x.train, y.train, 
                       type.measure="mse", alpha=1, family="gaussian")

alpha1.fit

#alpha1.predicted <- predict(alpha1.fit, s=alpha0.fit$lambda.1se, newx=x.test)

alpha1.predicted <- predict(alpha1.fit, s=alpha0.fit$lambda.min, newx=x.test)

#now calculate the mean square error of the true values verses the predicted values

alpha1.mse <- round(mean((y.test-alpha1.predicted)^2),digits=4)

alpha1.mse



# lets try something from snpStats

# determine which are the column indicies of the highest associated
# snps

#rs16849251
#rs13435610
#rs11466310
#rs13045127
#rs13345207

genotypes3_df
colnames(genotypes2)

colnames(genotypes[,20596])
which(colnames(genotypes2)=="rs13435610")


avg.mat3
imp <- snp.rhs.tests(genotypes3_mat~genotypes2[,20596], family="Gaussian", data=avg.mat3 )






# Run this once interactively to download and install BioConductor packages and other packages.


setwd("~/Dropbox/Ph.D/Bioinformatics/Project/datasets")


source ("http://bioconductor.org/biocLite.R ")
list.of.packages <- c("snpStats", "SNPRelate","rtracklayer", "biomaRt")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)


list.of.packages <- c('dplyr', 'ggplot2', 'coin' ,'igraph', 'devtools')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)



library(devtools)
install_url("http://cran.r-project.org/src/contrib/Archive/postgwas/postgwas_1.11.tar.gz")



source("globals.R")

# Downloading support files
# Download and unzip data needed for this tutorial


library(snpStats)

#*************************************************************
# load in the genotype data
file_prefix <- "ADNI_GO_2_Forward_Bin"
file.bed <- paste(getwd(),"/",file_prefix,".bed",sep="")
file.bim <- paste(getwd(),"/",file_prefix,".bim",sep="")
file.fam <- paste(getwd(),"/",file_prefix,".fam",sep="")

# read the data in plink format
geno <- read.plink(file.bed, file.bim, file.fam, na.strings = ("-9"))
str(ADNI.data)

#*************************************************************

# Read in PLINK files
#geno <- read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam, na.strings = ("-9"))

genotype <- geno$genotype
print(genotype)                  # 861473 SNPs read in for 1401 subjects



genoBim <- geno$map
colnames(genoBim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
print(head(genoBim))

rm(geno)



clinical <- read.csv("p_tau181_ADNI_GO_2.csv")

clinical.subcols <- clinical %>% select(RID, PLASMAPTAU181)

clinical.subcols

clinical

j <- 0
i <- 1

# avg.mat2 will contain the averages. rid.mat2 will contain the RID
avg.mat2 <-c()
rid.mat2 <-c()

while (i < 3758) {
  j <- 0
  while (pheno.subcols[i,'RID']==pheno.subcols[i+j,'RID']) {
    j = j+1
    print("equal")
    print(pheno.subcols[i+j,'RID'])
  }
  print("not equal")
  print(pheno.subcols[i+j,'RID'])
  
  average <- mean(pheno.subcols[(i:(i+j-1)),"PLASMAPTAU181"])
  rid <- pheno.subcols[(i+j-1),"RID"]
  i = i + j
  #print(average)
  #print(i)
  #print(j)
  avg.mat2 <- append(avg.mat2, average)
  rid.mat2 <- append(rid.mat2, rid)
}


averaged.clinical <- cbind(rid.mat2,avg.mat2)

# intersect the genotype data with the phenotype data
common_samples <- intersect(rownames(genotype),averaged.clinical[,1])


str(averaged.clinical)

# Subset genotype for subject data
genotype <- genotype[common_samples,]

print(genotype) # 98 subjects and 730525 SNPs


# Create SNP summary statistics (MAF, call rate, etc.)
snpsum.col <- col.summary(genotype)
print(head(snpsum.col))


# Setting thresholds
call <- 0.95
minor <- 0.01

# Filter on MAF and call rate
use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE                # Remove NA's as well

cat(ncol(genotype)-sum(use))
#    "54005 SNPs will be removed due to low MAF or call rate.

cat((ncol(genotype)-sum(use))/ncol(genotype))
#

ncol(genotype)
sum(use)


# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotype <- genotype[,use]
snpsum.col <- snpsum.col[use,]

print(genotype)   # 98 subjects and 676520 SNPs 


library(snpStats)
library(SNPRelate)               # LD pruning, relatedness, PCA
library(dplyr)

# Create sample statistics (Call rate, Heterozygosity)
snpsum.row <- row.summary(genotype)

# Add the F stat (inbreeding coefficient) to snpsum.row
MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(snpsum.row, Heterozygosity*(ncol(genotype))*Call.rate)
snpsum.row$hetF <- 1-(hetObs/hetExp)

head(snpsum.row)# 658186 SNPs remain


# Setting thresholds
sampcall <- 0.95    # Sample call rate cut-off
hetcutoff <- 0.1    # Inbreeding coefficient cut-off

sampleuse <- with(snpsum.row, !is.na(Call.rate) & Call.rate > sampcall & abs(hetF) <= hetcutoff)
sampleuse[is.na(sampleuse)] <- FALSE    # remove NA's as well
cat(nrow(genotype)-sum(sampleuse), 
    "subjects will be removed due to low sample call rate or inbreeding coefficient.\n") #0 subjects removed


# Subset genotype and clinical data for subjects who pass call rate and heterozygosity crtieria
genotype <- genotype[sampleuse,]
clinical<- clinical[ rownames(genotype), ]

# Checking for Relatedness

ld.thresh <- 0.2    # LD cut-off
kin.thresh <- 0.1   # Kinship cut-off

bed.fn <- "ADNI_GO_2_Forward_Bin.bed"
fam.fn <- "ADNI_GO_2_Forward_Bin.fam"
bim.fn <- "ADNI_GO_2_Forward_Bin.bim"

# Create gds file, required for SNPRelate functions
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "ADNIGO9.gds")


genofile <- openfn.gds("ADNIGO9.gds", readonly = FALSE)

# Automatically added "-1" sample suffixes are removed
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = TRUE)




#*************************[START: PRUNE IF SAMPLES ARE HIGHLY COORELATED]*************************

#Prune SNPs for IBD analysis
set.seed(1000)
geno.sample.ids <- rownames(genotype)

geno.sample.ids
show(index.gdsn(genofile, "sample.id"))
genofile$root[sample.id]

snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
                          #sample.id = geno.sample.ids, # Only analyze the filtered samples
                          snp.id = colnames(genotype)) # Only analyze the filtered SNPs
snpset.ibd <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.ibd),"will be used in IBD analysis\n")  #85225 will be used in IBD analysis

length(snpset.ibd)/ncol(genotype)

# Find IBD coefficients using Method of Moments procedure.  Include pairwise kinship.
ibd <- snpgdsIBDMoM(genofile, kinship=TRUE,
                    sample.id = geno.sample.ids,
                    snp.id = snpset.ibd,
                    num.thread = 1)


ibdcoeff <- snpgdsIBDSelection(ibd)     # Pairwise sample comparison
head(ibdcoeff)


# Check if there are any candidates for relatedness
ibdcoeff <- ibdcoeff[ ibdcoeff$kinship >= kin.thresh, ]

# iteratively remove samples with high kinship starting with the sample with the most pairings
related.samples <- NULL
while ( nrow(ibdcoeff) > 0 ) {
  
  # count the number of occurrences of each and take the top one
  sample.counts <- arrange(count(c(ibdcoeff$ID1, ibdcoeff$ID2)), -freq)
  rm.sample <- sample.counts[1, 'x']
  cat("Removing sample", as.character(rm.sample), 'too closely related to', 
      sample.counts[1, 'freq'],'other samples.\n')
  
  # remove from ibdcoeff and add to list
  ibdcoeff <- ibdcoeff[ibdcoeff$ID1 != rm.sample & ibdcoeff$ID2 != rm.sample,]
  related.samples <- c(as.character(rm.sample), related.samples)
}

# filter genotype and clinical to include only unrelated samples
genotype <- genotype[ !(rownames(genotype) %in% related.samples), ]
clinical <- clinical[ !(clinical$FamID %in% related.samples), ]

geno.sample.ids <- rownames(genotype)

cat(length(related.samples), 
    "similar samples removed due to correlation coefficient >=", kin.thresh,"\n") 


print(genotype)  

#*************************[END: PRUNE IF SAMPLES ARE HIGHLY COORELATED]*************************
# Checking for ancestry

geno.sample.ids

# Find PCA matrix
pca <- snpgdsPCA(genofile,
                 #sample.id = geno.sample.ids, 
                 snp.id = snpset.ibd, num.thread=1)


# Create data frame of first two principal comonents
pctab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],    # the first eigenvector
                    PC2 = pca$eigenvect[,2],    # the second eigenvector
                    stringsAsFactors = FALSE)

# Plot the first two principal comonents
plot(pctab$PC2, pctab$PC1, xlab="Principal Component 2", ylab="Principal Component 1", 
     main = "Ancestry Plot")


# Hardy-Weinberg SNP filtering on controls

hardy <- 10^-6      # HWE cut-off

CADcontrols
CADcontrols <- clinical[ clinical$CAD==0, 'FamID' ]
#snpsum.colCont <- col.summary( genotype[CADcontrols,] )
snpsum.colCont <- col.summary( genotype )
HWEuse <- with(snpsum.colCont, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
rm(snpsum.colCont)

HWEuse[is.na(HWEuse)] <- FALSE          # Remove NA's as well
cat(ncol(genotype)-sum(HWEuse),"SNPs will be removed due to high HWE.\n")  # 10313 SNPs will be removed due to high HWE

(ncol(genotype)-sum(HWEuse))/ncol(genotype)

# Subset genotype and SNP summary data for SNPs that pass HWE criteria
genotype <- genotype[,HWEuse]

print(genotype)              


# Save genotype and SNVs filtered data to use in later analyses
save(genotype, genoBim, clinical, file=working.data.fname(4))# 656890 SNPs remain



#Set LD threshold to 0.2
ld.thresh <- 0.2

set.seed(1000)
geno.sample.ids <- rownames(genotype)

genofile <- openfn.gds("ADNIGO9.gds", readonly = FALSE)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
                          sample.id = geno.sample.ids, # Only analyze the filtered samples
                          snp.id = colnames(genotype)) # Only analyze the filtered SNPs


snpset.pca <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.pca),"\n")  #76665  SNPs will be used in PCA analysis

pca <- snpgdsPCA(genofile,
                 #sample.id = geno.sample.ids,  
                 snp.id = snpset.pca, num.thread=1)


# Find and record first 10 principal components
# pcs will be a N:10 matrix.  Each column is a principal component.
pcs <- data.frame(FamID = pca$sample.id, pca$eigenvect[,1 : 10],
                  stringsAsFactors = FALSE)
colnames(pcs)[2:11]<-paste("pc", 1:10, sep = "")

print(head(pcs))

pcs$pc1



# Merge clincal data and principal components to create phenotype table
phenoSub <- merge(clinical,pcs)      # data.frame => [ FamID CAD sex age hdl pc1 pc2 ... pc10 ]


averaged.clinical[,2]

str(avg.mat2)

phenoSub$phenotype <- NA

norm_averaged.clinical <- log10(averaged.clinical[,2])

phenoSub$phenotype <- norm_averaged.clinical


# Show that the assumptions of normality met after transformation
par(mfrow=c(1,2))
hist(norm_averaged.clinical, main="Histogram of PTAU181", xlab="Frequency")


rownames(phenoSub) <- phenoSub$id

phenoSub$pc1

str(phenotype)
str(pcs$pc1)


pc1 <- pcs$pc1[common_samples]
pc2 <- pcs$pc2[common_samples]
pc3 <- pcs$pc3[common_samples]
pc4 <- pcs$pc4[common_samples]
pc5 <- pcs$pc5[common_samples]
pc6 <- pcs$pc6[common_samples]
pc7 <- pcs$pc7[common_samples]
pc8 <- pcs$pc8[common_samples]
pc9 <- pcs$pc9[common_samples]
pc10 <- pcs$pc10[common_samples]

str(phenotype)
str(pc1)


imp <- snp.rhs.tests(phenotype ~ pc1 + pc2 + pc3
                     + pc4 + pc5 + pc6
                     + pc7 + pc8 + pc9 + pc10,
                     family = "Gaussian", data = data.frame(phenotype) , snp.data = genotype)



results <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results <- results[!is.na(results$p.value),]

top.results <- results[(results$p.value)<1.5E-03,]

sort(top.results$p.value,descending=FALSE)

