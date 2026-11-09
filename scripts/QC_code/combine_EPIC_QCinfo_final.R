# Author: Eve Schoeffler (adapted from scripts from Yinan Zheng)


# Set libPaths and load libraries
.libPaths("path to R package library")
library(minfi)
library(purrr)

metric_file_path="path to directory containing metric files"
# Combine each metric

###### Bisulfite ############
print("Combining bisulfite data...")
# Set working directory
setwd(paste0(metric_file_path,"/QCsplitter/bisul/"))

# Set the number of plates
plates <- 28

# Create a list to hold the bisul chunks
bisul_list <- list.files(pattern = ".rds", recursive = T)
bisul_list <- unique(bisul_list)
length(bisul_list) #plates

bisul_data <- bisul_list %>% map(readRDS) # readRDS for each elem in bisul_list

bisul = NULL

# For each of the plates plates
for (i in 1:plates){ 
  bisul = c(bisul, bisul_data[[i]]) # Concatenating bisul vectors (for each sample on each plate)
}

rm(bisul_data)

# threshold of bisulfite conversion control intensity
bisulthre=mean(bisul,na.rm=TRUE)-3*sd(bisul,na.rm=TRUE)

# low qual samples based on bisulfite intensity
bisul_flag <- bisul < bisulthre
cat(sum(bisul_flag)," samples with bisulfite intensity less than ", bisulthre, "\n") # 19 samples

# Collect bad samples
bad_bisul_sample=names(bisul)[bisul_flag]

###### QC mat ################
print("Combining QC matrices...")
# Set working directory
setwd(paste0(metric_file_path,"/QCsplitter/qcmat/"))
# Create an empty list to hold each qcmatlist
qc_all <- vector(mode = "list", length = plates)

# Collect each plate's qcmatlist
for (plate in 1:plates) {
  # Set batch
  batch <- paste0("allEPIC_",plate,"_")
  
  # Take in all files
  qcmat_list <- list.files(pattern = ".rds", recursive = T)
  
  # Take in just the files from the current plate
  qcmat_list <- qcmat_list[grepl(paste0(batch), qcmat_list)]
  
  # Make sure the matlist is unique
  qcmat_list <- unique(qcmat_list)
  
  # Put the list in the proper order
  ordernum <- regmatches(gsub(paste0(batch), "", qcmat_list), gregexpr("[[:digit:]]+", gsub(paste0(batch), "", qcmat_list)))
  ordernum <- flatten_chr(ordernum)
  ordernum <-  as.numeric(ordernum)
  qcmat_list <- qcmat_list[order(ordernum)]
  #print(length(qcmat_list))
  qcmat_data <- qcmat_list %>% map(readRDS)
  qcmat_data <- as.matrix(do.call(rbind, qcmat_data))
  qc_all[[(plate)]] <- qcmat_data
}

# Transform to data.frame and create CpG column for merge
for (i in 1:length(qc_all)){
  qc_all[[i]] <- as.data.frame(qc_all[[i]])
  qc_all[[i]]$CpG <- rownames(qc_all[[i]])
}

# Check dimensions of each qc mat and check the missing
qcmat_dim <- vector(mode = "list", length = plates)
for (i in 1:length(qc_all)){
  qcmat_dim[[i]] <- c(dim(qc_all[[i]])[1], (dim(qc_all[[i]])[2]-1))
} 

# Merging all dataframes, only select intersect CpGs
qc_matrix <- Reduce(function(x,y) {merge(x,y, by = "CpG")}, qc_all) 
rownames(qc_matrix) <- qc_matrix$CpG
qc_matrix$CpG <- NULL
qc_matrix <- as.matrix(qc_matrix)
dim(qc_matrix) #865859   2709
qc_matrix[1:5,1:5]

###### Cross check bad samples ############
bad_samplethresample <- setdiff(names(bisul), colnames(qc_matrix)) 
length(bad_samplethresample) #7 (these are the samples that were removed based on samplethre)
setdiff(colnames(qc_matrix), names(bisul)) #character(0)

# # Check whether all the samples in bisul are also in qc_matrix
# common_samples <- intersect(names(bisul), colnames(qc_matrix))

# Remove the samples that failed samplethre from flag
flag <- bisul_flag[names(bisul_flag) %in% colnames(qc_matrix)]
length(flag) #2709 (was 2716)
cat(sum(flag)," samples with bisulfite intensity less than ", bisulthre, 
    " after removing samples that failed CpG threshold", "\n") # 16
# Remove samples that failed bisul from qc_matrix
qc_matrix <- qc_matrix[,!(colnames(qc_matrix) %in% bad_bisul_sample)]
dim(qc_matrix) # 865859   2693

badsample1=union(bad_bisul_sample, bad_samplethresample)
length(badsample1) #23 (7 samplethre + 19 bisul - 3 in common)
badsample1 <- badsample1[order(badsample1)]

# Confirm that the bad samples have been removed from qc_matrix
badsample2 <- setdiff(names(bisul), colnames(qc_matrix)) 
length(badsample2) #23
badsample2 <- badsample2[order(badsample2)]
identical(badsample1, badsample2) #TRUE

##### Define low quality CpGs ##############

sum(is.na(qc_matrix)) #0
NbadValuePerCpG <- apply(qc_matrix,1,sum)
length(NbadValuePerCpG) #865859
sum(is.na(NbadValuePerCpG)) #0

badValuePerCpG <- NbadValuePerCpG/ncol(qc_matrix)
CpGthre=0.05
flag2 <- badValuePerCpG>CpGthre & NbadValuePerCpG>1
cat(sum(flag2)," CpGs with percentage of low quality value greater than ",
    CpGthre,"\n") # 17839 

badCpG <- rownames(qc_matrix)[flag2]
length(badCpG) # 17839

qc_matrix=qc_matrix[!flag2,] 
dim(qc_matrix) # 848020   2693

rm(qc_all)

## I didn't look for bad CpGs based on dif annotations bc all our annotations
## should be the same

##### Identify outliers based on total intensity #############
print("Combining total intensity data...")
# Set working directory
setwd(paste0(metric_file_path,"/QCsplitter/mu/"))
# Create a list of all intensity chunks
mu_all <- vector(mode = "list", length = plates)

for (plate in 1:plates) {
  batch <- paste0("allEPIC_", plate, "_")
  mu_list <- list.files(pattern = ".rds", recursive = T)
  mu_list <- mu_list[grepl(paste0(batch), mu_list)]
  mu_all[[plate]] <- readRDS(mu_list)
}

# Only include the CpGs all batch have
mu_list <- vector(mode = "list", length = 4)
for (i in 1:length(mu_all)){
  mu_list[[i]] <- rownames(mu_all[[i]])
}

CpGs <- Reduce(intersect, mu_list)
#Reduce(function(x,y) {intersect(x,y)}, list)
length(CpGs) # 865859

for (i in 1:length(mu_all)){
  mu_all[[i]] <- mu_all[[i]][rownames(mu_all[[i]]) %in% CpGs,]
}

# Check again the dimension of each matrix
mu_dim <- vector(mode = "list", length = 4)
for (i in 1:length(mu_all)){
  mu_dim[[i]] <- dim(mu_all[[i]])
}
mu_dim

# Transform the mu list to a single matrix
mu_matrix <- as.matrix(do.call(cbind, mu_all))
dim(mu_matrix) # 865859   2709

# Remove bad CpGs and bad samples
dim(qc_matrix) # 848020   2693

mu_matrix1=mu_matrix[rownames(qc_matrix),]
mu_matrix1=mu_matrix1[,colnames(qc_matrix)]
dim(mu_matrix1) # 848020   2693

# Check mu matrix
# bad_bisul_sample <- badsample
bad_bisul_sample <- as.vector(t(bad_bisul_sample))

bad_samplethre_sample <- bad_samplethresample
bad_samplethre_sample <- as.vector(t(bad_samplethre_sample))
bad_sample <- c(bad_bisul_sample, bad_samplethre_sample) # not unique but shouldn't matter for operation

badCpG <- as.vector(t(badCpG))

mu_matrix2=mu_matrix[!(rownames(mu_matrix) %in% badCpG),]
mu_matrix2=mu_matrix2[,!(colnames(mu_matrix2) %in% bad_sample)]

dim(mu_matrix2) # 848020   2693

identical(dim(mu_matrix1), dim(mu_matrix2)) #TRUE

mu_matrix2 <- mu_matrix2[rownames(mu_matrix1),]
mu_matrix2 <- mu_matrix2[,colnames(mu_matrix1)]

identical(mu_matrix1, mu_matrix2) #TRUE

# Identify outliers
mumean=apply(mu_matrix1,2,mean,na.rm=TRUE)
length(mumean) # 3952
q2575 <- quantile(mumean, probs=c(0.25,0.75), na.rm=TRUE)
q2575
# 25%      75% 
# 8738.051 12216.131  
qr <- q2575["75%"]-q2575["25%"]
qr # 3478.079
low=q2575["25%"] - 3*qr
low # -1696.187
flag1=  mumean < low
cat("After excluding low quality samples and CpGs\n")
cat(sum(flag1)," samples are outliers based on averaged total intensity value","\n")
# 0 samples

##### Remove redundant objects to save memory ##################
ls()
rm(mu_all)
rm(mu_matrix)
rm(mu_matrix2)

mu <- mu_matrix1
rm(mu_matrix1)
gc()

##### Identify outliers based on beta values ###############
print("Combining beta values...")
# Set working directory
setwd(paste0(metric_file_path,"/QCsplitter/beta/"))
# Make list of all beta values
beta_all <- vector(mode = "list", length = plates)
for (plate in 1:plates) {
  batch <- paste0("allEPIC_", plate, "_")
  beta_list <- list.files(pattern = ".rds", recursive = T)
  beta_list <- beta_list[grepl(paste0(batch), beta_list)]
  beta_all[[plate]] <- readRDS(beta_list)
}

# Check dimension of each batch
beta_dim <- vector(mode = "list", length = plates)
for (i in 1:length(beta_all)){
  beta_dim[[i]] <- dim(beta_all[[i]])
}

identical(beta_dim, mu_dim) #TRUE

# Only include CpGs all batches have
list <- vector(mode = "list", length = plates)
for (i in 1:length(beta_all)){
  list[[i]] <- rownames(beta_all[[i]])
}

CpGs <- Reduce(intersect, list)
#Reduce(function(x,y) {intersect(x,y)}, list)
length(CpGs) # 865859

for (i in 1:length(beta_all)){
  beta_all[[i]] <- beta_all[[i]][rownames(beta_all[[i]]) %in% CpGs,]
}

# Transform beta list into single matrix
beta_matrix <- as.matrix(do.call(cbind, beta_all))
dim(beta_matrix) # 865859   2709

rm(beta_all)

# Remove bad CpGs and bad samples
beta_matrix=beta_matrix[rownames(qc_matrix),]
beta_matrix=beta_matrix[,colnames(qc_matrix)]
dim(beta_matrix) # 851282   3952

identical(rownames(qc_matrix), rownames(beta_matrix)) #TRUE
identical(colnames(qc_matrix), colnames(beta_matrix)) #TRUE

# Outliers in beta distribution
qq=apply(beta_matrix,2,function(x) quantile(x, probs=c(0.25,0.5,0.75), na.rm=TRUE))
q2575 <- apply(qq,1,function(x) quantile(x, probs=c(0.25,0.75), na.rm=TRUE))
q2575
#           25%       50%       75%
# 25% 0.14plates863 0.5868681 0.7669587
# 75% 0.2084576 0.7130019 0.8216196
qr <- q2575["75%",]-q2575["25%",]
qr
#       25%        50%        75% 
# 0.06417136 0.12613374 0.05466081
up=q2575["75%",] + 3*qr
up
#       25%       50%       75% 
# 0.4009717 1.0914031 0.9856020 
low=q2575["25%",] - 3*qr
low
#       25%        50%        75% 
# -0.04822783  0.20846689  0.60297631 
bflag=qq > up | qq < low
bflag=apply(bflag,2,sum)>0
cat(sum(bflag)," samples are outliers in beta value distribution","\n")
# 2  samples are outliers in beta value distribution
all_flag=bflag | flag1 # Flag if outlier in intensity or beta
outlier_sample=colnames(beta_matrix)[all_flag]
# [1] "206382230116_R01C01" "205841320073_R06C01"

cat(sum(all_flag)," outlier samples were added into badsample list\n")
#2  outlier samples were added into badsample list

outlier_sample <- as.vector(t(outlier_sample))
###### Reassemble QC ##########################
detPthre=0.000001
detPtype="negative"
nbthre=3
samplethre=0.05
CpGthre=0.05

length(bad_bisul_sample) #19
length(bad_samplethre_sample) #7
length(outlier_sample) #3

badsample <- c(bad_bisul_sample, bad_samplethre_sample, outlier_sample)
length(badsample) #29
sum(duplicated(badsample)) #3

badsample <- unique(badsample)
length(badsample) #26

badCpG <- as.vector(t(badCpG))
length(badCpG) # 17839

# qc = list(badsample=badsample,badCpG=badCpG,outlier_sample=outlier_sample)

##### Read in nbead and create a list ###################
print("Combining nbead matrices...")
# Set working directory
setwd(paste0(metric_file_path,"/QCsplitter/nbead/"))
nbead_all <- vector(mode = "list", length = plates)
for (plate in 1:plates) {
  batch <- paste0("allEPIC_", plate,"_")
  nbead_list <- list.files(pattern = ".rds", recursive = T)
  nbead_list <- nbead_list[grepl(paste0(batch), nbead_list)]
  nbead_all[[(plate)]] <- readRDS(nbead_list)
}

# Check dimension of each batch
nbead_dim <- vector(mode = "list", length = plates)
for (i in 1:length(nbead_all)){
  nbead_dim[[i]] <- dim(nbead_all[[i]])
}

# Only include the CpGs all batch have
list <- vector(mode = "list", length = plates)
for (i in 1:length(nbead_all)){
  list[[i]] <- rownames(nbead_all[[i]])
}

CpGs <- Reduce(intersect, list)
length(CpGs) # 865859

for (i in 1:length(nbead_all)){
  nbead_all[[i]] <- nbead_all[[i]][rownames(nbead_all[[i]]) %in% CpGs,]
}

# Transform nbead list into a single matrix
nbead <- as.matrix(do.call(cbind, nbead_all))
dim(nbead) #865859   2716

rm(nbead_all)
gc()

###### Read in detP and create a list
print("Combining detP matrices...")
# Set working directory
setwd(paste0(metric_file_path,"/QCsplitter/detP/"))
detP_all <- vector(mode = "list", length = plates)
for (plate in 1:plates) {
  batch <- paste0("allEPIC_", plate, "_")
  detP_list <- list.files(pattern = ".rds", recursive = T)
  detP_list <- detP_list[grepl(paste0(batch), detP_list)]
  detP_all[[plate]] <- readRDS(detP_list)
}

# Check dimension of each batch
detP_dim <- vector(mode = "list", length = plates)
for (i in 1:length(detP_all)){
  detP_dim[[i]] <- dim(detP_all[[i]])
}

identical(nbead_dim, detP_dim) # TRUE

# Only include CpGs all plates have
list <- vector(mode = "list", length = plates)
for (i in 1:length(detP_all)){
  list[[i]] <- rownames(detP_all[[i]])
}

CpGs <- Reduce(intersect, list)
#Reduce(function(x,y) {intersect(x,y)}, list)
length(CpGs) #865859

for (i in 1:length(detP_all)){
  detP_all[[i]] <- detP_all[[i]][rownames(detP_all[[i]]) %in% CpGs,]
}

# Transform detP list into single matrix
detP <- as.matrix(do.call(cbind, detP_all))
dim(detP) #865859   3966

rm(detP_all)

###### Combine all QC components ##################
qc = list(detP=detP,nbead=nbead,bisul=bisul,badsample=badsample,badCpG=badCpG,outlier_sample=outlier_sample)
rm(detP)
rm(nbead)
rm(bisul)
rm(badsample)
rm(badCpG)
rm(outlier_sample)
print("Saving final QC info...")
dir.create(paste0(metric_file_path,"/Final_QC"))
saveRDS(qc, paste0(metric_file_path,"/Final_QC/qc_info.EPIC.saliva.RDS"))
print("Done!")
