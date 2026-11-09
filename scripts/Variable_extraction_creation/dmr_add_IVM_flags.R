dmr_filter_fun<-function(dmrfilepath,ewasfilepath,samplecount){
  # flag files
  sex<-read.csv("path to Sex_DMP_Y9.csv", skip = 1)
  sex<-sex$Probe_ID
  man_epic <- readRDS("path to EPIC.hg19.manifest.dataframe.rds")
  epic_mask_general <- rownames(man_epic)[man_epic$MASK_general == T]
  tiss <- readRDS("path to saliva_blood_difMethCpGs.RDS")
  gap_probes<-readRDS("path to nmode.saliva.030926.RDS")
  gap_probes<-gap_probes[gap_probes$nmode>1,"probe"]
  blood_gap_probes<-readRDS("path to nmode.blood.021425.RDS")
  blood_gap_probes<-names(blood_gap_probes[blood_gap_probes>1])
  if(grepl("y22",dmrfilepath==TRUE)){
    flagged_probes<-list(Probes_mask_general=epic_mask_general,Probes_Sex_Difference=sex,Probes_nmode=blood_gap_probes,Probes_Saliva_Blood_Differencee=tiss)} else {
      flagged_probes<-list(Probes_mask_general=epic_mask_general,Probes_Sex_Difference=sex,Probes_nmode=gap_probes,Probes_Saliva_Blood_Differencee=tiss)} 
  
  ewas <- readRDS(ewasfilepath)
  
 if(!"Name" %in% names(ewas)){
                    ewas$Name<-ewas$ID}
  if(!"probe" %in% names(ewas)){
    ewas$probe<-ewas$ID}
  
  if(!"bacon.se_fixed" %in% names(ewas)){
    ewas$bacon.se_fixed<-ewas$bacon.se}
  if(!"bacon.est_fixed" %in% names(ewas)){
    ewas$bacon.est_fixed<-ewas$bacon.est}
  ewas$n_region_input <- samplecount$n[
    match(ewas$Name, samplecount$probe)
  ]
  ewas <- ewas[order(ewas$seqnames, ewas$start), ]
  
 
  
  
  anno_df<-readRDS(dmrfilepath)
  

calculate_dmr_ivm <- function(chr, start, end, ewas) {
  
  idx <- which(
    as.character(ewas$seqnames) == chr &
      ewas$start >= start &
      ewas$start <= end
  )
  
  se <- ewas$bacon.se_fixed[idx]
  est <- ewas$bacon.est_fixed[idx]
  n <- ewas$n_region_input[idx]
  
  valid <- is.finite(se) &
    is.finite(est) &
    is.finite(n) &
    se > 0 &
    n > 0
  
  if (!any(valid)) {
    return(c(
      IVM_coef = NA_real_,
      IVM_se = NA_real_,
      Nregion = NA_real_,
      n_CpGs = 0
    ))
  }
  
  w <- 1 / se[valid]^2
  
  c(
    IVM_coef = weighted.mean(est[valid], w),
    IVM_se = sqrt(1 / sum(w)),
    Nregion = weighted.mean(n[valid], w),
    n_CpGs = sum(valid)
  )
}

ivm_results <- t(vapply(
  seq_len(nrow(anno_df)),
  function(i) {
    calculate_dmr_ivm(
      chr = as.character(anno_df$seqnames[i]),
      start = anno_df$start[i],
      end = anno_df$end[i],
      ewas = ewas
    )
  },
  numeric(4)
))
anno_df <- cbind(anno_df, ivm_results)
anno_df$t_region<-anno_df$IVM_coef/anno_df$IVM_se
anno_df$df_region<-anno_df$Nregion-9-1
anno_df$partialR2<-(anno_df$t_region^2)/(anno_df$t_region^2+anno_df$df_region)



for(flag in names(flagged_probes)) {
  anno_df[[paste0("n_",flag)]] <-vapply(
    seq_len(nrow(anno_df)),
    function(i) {
      idx <-(
        as.character(ewas$seqnames) ==
          as.character(anno_df$seqnames[i]) &
        
          (ewas$start >= anno_df$start[i]) &
          (ewas$start <=anno_df$end[i])
      )
      probes <-unique(as.character(ewas$Name[idx]))
      sum(probes %in% flagged_probes[[flag]])
    },
    integer(1)
  )
}



return(anno_df)
}


# Age 9
y9metalist<-readRDS("path to y9_meta_analysis_odel1_ewas_output.RDS")

y9epicewas<-read.delim("path to y9_epic_model1.ewasoutput.txt")
y9epicewas<-y9epicewas[y9epicewas$probe %in% names(y9metalist),]
y9450Kewas<-read.delim("path to y9_450K_smodel1.ewasoutput.txt")

y9450Kewas<-y9450Kewas[y9450Kewas$probe %in% names(y9metalist),]
y9450Kewaso<-y9450Kewas[match(y9epicewas$probe,y9450Kewas$probe),]
y9samplecount<-data.frame(probe=y9450Kewaso$probe,hm450=as.numeric(y9450Kewaso$n),epic=as.numeric(y9epicewas$n))
y9samplecount$n<-y9samplecount$hm450+y9samplecount$epic

age9<-dmr_filter_fun(dmrfilepath = "path to y9_DMRs_DMRcate_combp_chipseeker.RDS",
               ewasfilepath="path to annotated y9 meta-analysis model 1 ewas output",samplecount = y9samplecount)


# Age 15
y15metalist<-readRDS("path to y15_meta_analysis_odel1_ewas_output.RDS")

y15epicewas<-read.delim("path to y15_epic_model1.ewasoutput.txt")
y15epicewas<-y15epicewas[y15epicewas$probe %in% names(y15metalist),]
y15450Kewas<-read.delim("path to y15_450K_model1.ewasoutput.txt")

y15450Kewas<-y15450Kewas[y15450Kewas$probe %in% names(y15metalist),]
y15450Kewaso<-y15450Kewas[match(y15epicewas$probe,y15450Kewas$probe),]
y15samplecount<-data.frame(probe=y15450Kewaso$probe,hm450=as.numeric(y15450Kewaso$n),epic=as.numeric(y15epicewas$n))
y15samplecount$n<-y15samplecount$hm450+y15samplecount$epic

age15<-dmr_filter_fun(dmrfilepath = "path to y15_DMRs_DMRcate_combp_chipseeker.RDS",
                         samplecount=y15samplecount,
                         ewasfilepath="path to y15_meta_analysis_odel1_ewas_output.RDS")
age15_filt<-age15[age15$dmrcate_present%in% "TRUE" & age15$combp_present %in% TRUE,]

# Mixed effect model - Age 9/Age15
y9y15dmr<-readRDS("path to y9y15_DMRs_DMRcate_combp_chipseeker.RDS")
y9y15ewas<-readRDS("path to y9_y15_ewas_fixed_effects.RDS")
y9y15ewassamplecount<-data.frame(probe=y9y15ewas$ID,hm450=as.numeric(y9y15ewas$n_obs_450k),epic=as.numeric(y9y15ewas$n_obs_epic))
y9y15ewassamplecount$n<-y9y15ewassamplecount$hm450+y9y15ewassamplecount$epic
y9y15<-dmr_filter_fun(dmrfilepath = "path to y9y15_DMRs_DMRcate_combp_chipseeker.RDS",ewasfilepath = 
                        "path to y9_y15_ewas_fixed_effects.RDS",samplecount = y9y15ewassamplecount)

# Age 22

y22<-readRDS("path to y22_epic_model_1_ewas_output.RDS")
y22samplecount<-y22[,c("probe","n")]

age22<-dmr_filter_fun(ewasfilepath = "path to y22_epic_model_1_ewas_output.RDS",
               dmrfilepath="path to 22_DMRs_DMRcate_combp_chipseeker.RDS",samplecount = y22samplecount)
#### Save files
write.table(age9,file="age9_DMRs_annotated_IVM.txt",sep="\t",quote=F,row.names=F)
write.table(age15,file="age15_DMRs_annotated_IVM.txt",sep="\t",quote=F,row.names=F)
write.table(age22,file="age22_DMRs_annotated_IVM.txt",sep="\t",quote=F,row.names=F)
write.table(y9y15,file="age9age15_DMRs_annotated_IVM.txt",sep="\t",quote=F,row.names=F)
