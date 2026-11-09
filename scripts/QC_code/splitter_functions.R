# Author: Eve Schoeffler (adapted from scripts from Yinan Zheng)

function_list <- c("matsplitter", "QCsplitter")

trash <- sapply(function_list, function(x) suppressWarnings(rm(x)))

message("The following tools have been succesfully loaded:")
trash <- sapply(function_list, function(x) message(x))

##Function to split matrix into chunks###
matsplitter <- function(mat, n, direction = "rows") { # Give a matrix, number of chunks, and a direction
  
  # If splitting the matrix by rows (basically make into smaller rgsets)
  if (direction == "rows") {
    
    # Calculate the number of rows in each chunk
    rg <- (dim(mat)[1])%/%n   # %/% ensures the output of division is an integer
    rg1 <- (dim(mat)[1])%/%n + (dim(mat)[1])%%n # Number of rows for last chunk (+ %% n so that it includes remainder)
    #(n-1)*rg + rg1 == dim(mat)[1]
    
    # Create a list to hold matrices for each chunk
    my_list <- vector(mode = "list", length = n)
    for (i in 1:(n-1)){ # Head of for-loop
      new_element <- mat[((i-1)*rg+1):(i*rg), ]    # Select rows for each sub-matrix
      my_list[[i]] <- new_element    # Fill each list element
    } 
    
    # Create the final element in the list
    my_list[[n]] <- mat[((n-1)*rg+1):(dim(mat)[1]), ]
    message("Exported ", length(my_list), " chunks of matrix, split across ", direction)
    
  # If splitting the matrix by columns (splitting different CpGs apart, keeping all samples)
  } else if (direction == "cols") {
    
    # Calculate the number of columns in each chunk
    cg <- (dim(mat)[2])%/%n
    cg1 <- (dim(mat)[2])%/%n + (dim(mat)[2])%%n
    #(n-1)*cg + cg1 == dim(mat)[2]
    
    # Create a list to hold matrices for each chunk
    my_list <- vector(mode = "list", length = n)
    for (i in 1:(n-1)) { # Head of for-loop
      new_element <- mat[ ,((i-1)*cg+1):(i*cg)]        # Select columns for each sub-matrix
      my_list[[i]] <- new_element    # Fill each list element
    } 
    
    # Create the final element in the list
    my_list[[n]] <- mat[ ,((n-1)*cg+1):(dim(mat)[2])]
    message("Exported ", length(my_list), " chunks of matrix, split across ", direction)
    
  } else {stop("direction should be 'rows' or 'cols'!")}
  
  # Return the list of matrices
  return(my_list)
}

# ###Testing of the function####
# mylist <- matsplitter(qcmat, 20, direction = "rows")
# alldata <- as.matrix(do.call(rbind, mylist))

# mylist <- matsplitter(qcmat, 5, direction = "cols")
# alldata <- as.matrix(do.call(cbind, mylist))

# ###Test if dims are identical for recombined matrix and original one
# dim(qcmat) #866836     11
# dim(alldata) #866836     11
# identical(dim(qcmat), dim(alldata)) #TRUE
# ###Test if CpGs are identical for recombined matrix and original one
# identical(row.names(qcmat), row.names(alldata))  #TRUE
# ###Test if two matrices are identical
# identical(qcmat, alldata) #TRUE


###Function to read in rgDataSet or RGChannelSetExtended and output chunks of QC elements
QCsplitter <- function(rgSet,detPthre=0.000001,detPtype="negative",nbthre=3,samplethre=0.05,CpGthre=0.05,
                       bisulthre=NULL,outlier=TRUE, chunk_num=10)
{
  
  ##Load packages
  library(minfi)
  library(ENmix)
  ##number of bead
  if(!is(rgSet, "rgDataSet") & !is(rgSet, "RGChannelSetExtended"))
    stop("[QCinfo] The input should be an object of 'rgDataSet' or 'RGChannelSetExtended'")
  
  if(is(rgSet, "rgDataSet")){
    cginfo=getCGinfo(rgSet)
    typeI <- cginfo[cginfo$Infinium_Design_Type %in% c("I","snpI"),]
    typeIred=typeI[typeI$Color_Channel=="Red",]
    typeIgrn=typeI[typeI$Color_Channel=="Grn",]
    typeII<-cginfo[cginfo$Infinium_Design_Type %in% c("II","snpII"),]
    locusNames=c(typeIred$Name,typeIgrn$Name,typeII$Name)
    ##detection P value
    detP<-calcdetP(rgSet,detPtype=detPtype)
    ctrls<-metadata(rgSet)$ictrl
    
  # Check if the object is an RGChannelSetExtended and access information  
  }else if(is(rgSet, "RGChannelSetExtended")){
    
    typeI <- getProbeInfo(rgSet, type = "I")
    typeII <- getProbeInfo(rgSet, type = "II")
    locusNames <- getManifestInfo(rgSet, "locusNames")
    
    ##detection P value
    detP<-detectionP(rgSet)
    ctrls<-getProbeInfo(rgSet,type="Control")
  }
  
  # Get the number of Type I beads at AddressA
  bc_I <- assays(rgSet)$NBeads[typeI$AddressA,]
  
  # Flag if there are more beads at AddressA than AddressB
  flag<-bc_I>assays(rgSet)$NBeads[typeI$AddressB,]
  
  # Change the number of Type I beads to the smaller option (between A and B)
  bc_I[flag] <- assays(rgSet)$NBeads[typeI$AddressB,][flag];
  
  # Create an empty nbead matrix which has all samples in rgSet and all loci
  nbead <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
                  dimnames = list(locusNames, colnames(rgSet)))
  
  # At the rows that are type I probes, take the nBeads values from bc_I
  nbead[typeI$Name,]<-bc_I
  
  # At the rows that are type II probes, get the number of beads at Address A (I guess it's consistent for type II)
  nbead[typeII$Name,]<-assays(rgSet)$NBeads[typeII$AddressA,]
  rm(list=c("bc_I","flag"))
  
  # Get rid of extra rows and columns that do not overlap between detP matrix and nbead matrix
  if(!identical(rownames(detP),rownames(nbead))){
    idx=intersect(rownames(detP),rownames(nbead))
    detP=detP[idx,]
    nbead=nbead[idx,]
  }
  if(!identical(colnames(detP),colnames(nbead))){
    idx=intersect(colnames(detP),colnames(nbead))
    detP=detP[,idx]
    nbead=nbead[,idx]
  }
  
  ## bisulfite conversion internal controls
  
  # Take only the controls which are included in rgSet
  ctrls=ctrls[ctrls$Address %in% rownames(rgSet),]
  ctrl_r <- assays(rgSet)$Red[ctrls$Address,] # Extract control red
  ctrl_g <- assays(rgSet)$Green[ctrls$Address,] # Extract control green
  cc=ctrls[(ctrls$Type %in% c("BISULFITE CONVERSION I")) & ((ctrls$Color %in% 
                                                               c("Green","Lime","LimeGreen")) | (ctrls$ExtendedType %in% 
                                                                                                   c("ctl-BISULFITE-CONVERSION-I-140M_MUS","ctl-BISULFITE-CONVERSION-I-303M_MUS"))),]
  I_green=colMeans(ctrl_g[cc$Address,])
  cc=ctrls[(ctrls$Type %in% c("BISULFITE CONVERSION I")) & ((ctrls$Color %in%
                                                               c("Purple","Red","Tomato")) | (ctrls$ExtendedType %in% 
                                                                                                c("ctl-BISULFITE-CONVERSION-I-318M_MUS","ctl-BISULFITE-CONVERSION-I-330U_MUS"))),]
  I_red=colMeans(ctrl_r[cc$Address,])
  cc=ctrls[ctrls$Type %in% c("BISULFITE CONVERSION II") & ctrls$Color %in% c("Crimson",
                                                                             "DarkMagenta","Red","Orange","Purple","Tomato"),]
  II_red=colMeans(ctrl_r[cc$Address,])
  bisul=(I_green+I_red+II_red)/3
  
  # #threshold of bisulfite conversion control intensity
  # if(is.null(bisulthre)){bisulthre=mean(bisul,na.rm=TRUE)-3*sd(bisul,na.rm=TRUE)}
  
  ##low quality samples
  
  # Create a logical matrix to show if CpGs fail QC
  qcmat <- nbead<nbthre | detP>detPthre
  
  # Calculate proportion of failed QC CpGs per samples (column)
  badValuePerSample <- apply(qcmat,2,sum)/nrow(qcmat)
  
  # Flag samples that fall below threshold
  flag <- badValuePerSample > samplethre ##| bisul < bisulthre##
  message(sum(flag)," samples with percentage of low quality CpG value greater than ",
          samplethre)
  
  # Create a list of bad samples
  badsample=colnames(qcmat)[flag]
  
  # Keep only high quality samples
  qcmat <- qcmat[,!flag]
  
  ##split qc matrix into chunk_num chunks across rows (CpGs)
  qcmatlist <- matsplitter(qcmat, chunk_num, direction = "rows")
  
  if(outlier){
    if(is(rgSet, "rgDataSet")){mdat=getmeth(rgSet)
    }else if(is(rgSet, "RGChannelSetExtended")){mdat=preprocessRaw(rgSet)}
    rm(rgSet)}
  
  #Identifying outlier samples
  if(outlier){
    message("Identifying ourlier samples based on beta or total intensity values...")
    mdat=mdat[rownames(qcmat),]
    mdat=mdat[,colnames(qcmat)]
    #outliers based on total intensity values
    mu <- assays(mdat)$Meth+assays(mdat)$Unmeth
    #outliers in beta value distribution
    beta=getB(mdat, type="Illumina") 
  }
  
  rm(mdat)
  
  if(outlier)
  {list(detP=detP,nbead=nbead,bisul=bisul, badsample=badsample, qcmatlist=qcmatlist,
        mu=mu, beta=beta)
  }else{list(detP=detP,nbead=nbead,bisul=bisul, badsample=badsample, qcmatlist=qcmatlist)}
}

# ###Testing
# qc <- QCsplitter(rgSet,detPthre=0.000001,detPtype="negative",nbthre=3,samplethre=0.05,CpGthre=0.05,
#                        bisulthre=NULL,outlier=TRUE, chunk_num=10)