library(GenomicRanges)

setwd("path to working directory")

OUTDIR="/drives/drive1/ff_chaya/Lisa/April_2026/Sept2026_revision"
dir.create<-OUTDIR


y9dmrff_file<-"path to dmrff_y9_sig_meta_output_chipseeker.RDS"
y9dmrcate_file<-"path to y9_meta_dmrcate_chipseeker.RDS"
y9combp_file<-"path to y9_meta_combp_regions.1500-200_chipseeker_significanthits.txt"

y15dmrff_file<-"path to dmrff_y15_sig_meta_output_chipseeker.RDS"
y15dmrcate_file<-"path to y15_dmrcate_chipseeker.RDS"
y15combp_file<-"path to y15_meta_combp_regions.1500-200_chipseeker_significanthits.txt"

y22dmrff_file<-"path to y22_epic_dmrff_ewasoutput_sig_chipseeker.RDS"
y22dmrcate_file<-"path to y22_epic_dmrcate_chipseeker.RDS"
y22combp_file<-"path to y22_epic_combpinput.regions.sig.1500-200_chipseeker_significanthits_partialR2.txt"

y9_y15dmrff_file<-"path to y9_y15_epic_dmrff_ewasoutput_sigdmrffoutput.txt"
y9_y15dmrcate_file<-"path to y9_15_longitudinal_blood_dmrcate_resultsdf.RDS"
y9_y15combp_file<-"path to y9_y15_ewas_fixed_effects.baconadj.combpinput.regions.sig.1500-200_chipseeker_significanthits_partialR2.txt"
#

compare_dmrs_by_method_fun<-function(dmrff_file, dmrcate_file, combp_file){
  dmrff<-if(grepl("rds",dmrff_file,ignore.case=T)==TRUE) {readRDS(dmrff_file)} else {read.delim(dmrff_file)}
  #names(dmrff)<-gsub("\\.","_",names(dmrff))
  if(!"dmrff_chr" %in% names(dmrff)){names(dmrff)[names(dmrff)%in% "dmr.chr"]<-"dmrff_chr"} else {print("column names have dmrff prefix")}
  if(!"dmrff_start" %in% names(dmrff)){names(dmrff)[names(dmrff)%in% "dmr.start"]<-"dmrff_start"} else {print("column names have dmrff prefix")}
  if(!"dmrff_end" %in% names(dmrff)){names(dmrff)[names(dmrff)%in% "dmr.end"]<-"dmrff_end"} else {print("column names have dmrff prefix")}
  if(!"dmrff_estimate" %in% names(dmrff)){names(dmrff)[names(dmrff)%in% "dmr.estimate"]<-"dmrff_estimate"} else {print("column names have dmrff prefix")}
  if(!"chr" %in% dmrff$dmrff_chr){dmrff$dmrff_chr<-paste0("chr",dmrff$dmrff_chr)}else {print("dmrff_chr already has chr prefix")}
  if(!"dmrff_estimate" %in% names(dmrff)){names(dmrff)[names(dmrff) =="estimate"]<-c("dmrff_estimate")}else {print("column names have dmrff prefix")}

  dmrff<-dmrff[,c("dmrff_chr","dmrff_start","dmrff_end","dmrff_estimate")]
  dmrff<-dmrff[!duplicated(dmrff),]
  names(dmrff)[4]<-"estimate"
  
  dmrcate<-if(grepl("rds",dmrcate_file,ignore.case=T)==TRUE) {readRDS(dmrcate_file)} else {read.delim(dmrcate_file)}
  dmrcate$estimate<-dmrcate$meandiff
  dmrcate<-dmrcate[,c("seqnames","start","end","estimate")]
  names(dmrcate)[4]<-"estimate"
  
  combp<-if(grepl("rds",combp_file,ignore.case=T)==TRUE) {readRDS(combp_file)} else {read.delim(combp_file)}
  combp<-combp[,c("seqnames","start","end","IVM_coef")]
  names(combp)[4]<-"estimate"
  
  dmrlist<-list()
  dmrlist$dmrff<-makeGRangesFromDataFrame(dmrff,seqnames="dmrff_chr",start.field="dmrff_start",end.field = "dmrff_end",keep.extra.columns = T)
  dmrlist$dmrcate<-makeGRangesFromDataFrame(dmrcate,seqnames="seqnames",start.field="start",end.field = "end",keep.extra.columns = T)
  dmrlist$combp<-makeGRangesFromDataFrame(combp,seqnames="seqnames",start.field="start",end.field = "end",keep.extra.columns = T)
  return(dmrlist)}

make_dmr_master <- function(dmr_list,
                              estimate_col = "estimate") {
    
    comparison_names<-names(dmr_list)
    gr1<-dmr_list[[1]]
    gr2<-dmr_list[[2]]
    gr3<-dmr_list[[3]]
    
    
    ## 1. Create master set of collapsed genomic regions
    master <- reduce(c(gr1, gr2, gr3))
    master$region_id <- paste0("DMR_", seq_along(master))
    
    ## 2. Add information from each original DMR set
    for (i in seq_along(dmr_list)) {
      dmr <- dmr_list[[i]]
      prefix <- comparison_names[i]
      ## Find overlaps between collapsed regions and original DMRs
      hits <- findOverlaps(master, dmr)
      
      
      ## Presence / absence
      present <- logical(length(master))
      present[unique(queryHits(hits))] <- TRUE
      mcols(master)[[paste0(prefix, "_present")]] <- present
      
      ## Original coordinates
     coords <- paste0(as.character(seqnames(dmr)), ":",start(dmr), "-",end(dmr))
      
      coord_by_region <- split(coords[subjectHits(hits)],queryHits(hits))
      original_coords <- rep(NA_character_, length(master))
      original_coords[as.integer(names(coord_by_region))] <-vapply(coord_by_region,
          function(x)paste(unique(x), collapse = "; "),character(1))
      mcols(master)[[paste0(prefix, "_coordinates")]] <-original_coords

      ## Original estimate(s)
      est <- mcols(dmr)[[estimate_col]]
      est_by_region <- split(est[subjectHits(hits)],queryHits(hits))
      estimates <- rep(NA_character_, length(master))
      estimates[as.integer(names(est_by_region))] <-vapply(est_by_region,function(x) paste(round(x, 4), collapse = "; "),
          character(1))
      mcols(master)[[paste0(prefix, "_estimate")]] <-estimates
    
    ## Direction
    direction_by_region <- lapply(
      est_by_region,
      function(x) {
        
        dirs <- ifelse(
          x > 0, "Hyper",
          ifelse(x < 0, "Hypo", "No_change")
        )
        
        unique(dirs)
      }
    )
    directions <- rep(NA_character_, length(master))
    directions[as.integer(names(direction_by_region))] <-
      vapply(direction_by_region,function(x)paste(x, collapse = "; "),character(1))
    mcols(master)[[paste0(prefix, "_direction")]] <-directions
  
    ## Coordinate + estimate + direction to keep everything paired
    original_direction <- ifelse(est > 0, "Hyper",ifelse(est < 0, "Hypo", "No_change"))
    detail <- paste0(coords," [estimate=",round(est, 4),", ",original_direction,"]")
    detail_by_region <- split(detail[subjectHits(hits)],queryHits(hits))
    details <- rep(NA_character_, length(master))
    details[as.integer(names(detail_by_region))] <-vapply(detail_by_region,function(x)
          paste(x, collapse = "; "),character(1))
        mcols(master)[[paste0(prefix, "_details")]] <-details
  }
    # 3. Assign Venn category
    present_cols <- paste0(comparison_names, "_present")
    presence <- as.data.frame(mcols(master)[, present_cols, drop = FALSE])
    master$category <- apply(presence,1,function(x) {paste(comparison_names[x],collapse = "_")})
    
    ## 4. Number of DMR sets containing each region
    master$n_comparisons <- rowSums(presence)
    
    ## Return master GRanges
    return(master)
  }

y9<-compare_dmrs_by_method_fun(dmrff_file = y9dmrff_file,dmrcate_file = y9dmrcate_file,combp_file = y9combp_file)
y9_regions <- make_dmr_master(dmr_list = y9)
y9_regionsdf <- as.data.frame(y9_regions)
y9_any2<-y9_regionsdf[y9_regionsdf$n_comparisons>=2,]
saveRDS(y9_any2,file=paste0(OUTDIR,"/y9_DMRs_any2methods.RDS"))
y9_dmrcate_combp<-y9_regionsdf[y9_regionsdf$dmrcate_present==TRUE & y9_regionsdf$combp_present==TRUE,]
saveRDS(y9_dmrcate_combp,file=paste0(OUTDIR,"/y9_DMRs_DMRcate_combp.RDS"))

y15<-compare_dmrs_by_method_fun(dmrff_file = y15dmrff_file,dmrcate_file = y15dmrcate_file,combp_file = y15combp_file)
y15_regions <- make_dmr_master(dmr_list = y15)
y15_regionsdf <- as.data.frame(y15_regions)
y15_any2<-y15_regionsdf[y15_regionsdf$n_comparisons>=2,]
y15_dmrcate_combp<-y15_regionsdf[y15_regionsdf$dmrcate_present==TRUE & y15_regionsdf$combp_present==TRUE,]
saveRDS(y15_any2,file=paste0(OUTDIR,"/y15_DMRs_any2methods.RDS"))
saveRDS(y15_dmrcate_combp,file=paste0(OUTDIR,"/y15_DMRs_DMRcate_combp.RDS"))

y22<-compare_dmrs_by_method_fun(dmrff_file = y22dmrff_file,dmrcate_file = y22dmrcate_file,combp_file = y22combp_file)
y22_regions <- make_dmr_master(dmr_list = y22)
y22_regionsdf <- as.data.frame(y22_regions)
y22_any2<-y22_regionsdf[y22_regionsdf$n_comparisons>=2,]
y22_dmrcate_combp<-y22_regionsdf[y22_regionsdf$dmrcate_present==TRUE & y22_regionsdf$combp_present==TRUE,]
saveRDS(y22_any2,file=paste0(OUTDIR,"/y22_DMRs_any2methods.RDS"))
saveRDS(y22_dmrcate_combp,file=paste0(OUTDIR,"/y22_DMRs_DMRcate_combp.RDS"))

y9y15<-compare_dmrs_by_method_fun(dmrff_file = y9_y15dmrff_file,dmrcate_file = y9_y15dmrcate_file,combp_file = y9_y15combp_file)
y9y15_regions <- make_dmr_master(dmr_list = y9y15)
y9y15_regionsdf <- as.data.frame(y9y15_regions)
y9y15_any2<-y9y15_regionsdf[y9y15_regionsdf$n_comparisons>=2,]
saveRDS(y9y15_any2,file=paste0(OUTDIR,"/y9y15_DMRs_any2methods.RDS"))
y9y15_dmrcate_combp<-y9y15_regionsdf[y9y15_regionsdf$dmrcate_present==TRUE & y9y15_regionsdf$combp_present==TRUE,]
saveRDS(y9y15_dmrcate_combp,file=paste0(OUTDIR,"/y9y15_DMRs_DMRcate_combp.RDS"))


persistentdmr_function2<-function(y9,y15,y22){
  y9gr<-makeGRangesFromDataFrame(y9)
  y15gr<-makeGRangesFromDataFrame(y15)
  y22gr<-makeGRangesFromDataFrame(y22)
  y9y15gr<-intersect(y9gr,y15gr)
  y9y15y22<-intersect(y9y15gr,y22gr)
  return(as.data.frame(y9y15y22))
}

persistent_dmrs_fun<-
function(y9,y15,y22){
  # Venn diagram
  y9gr<-makeGRangesFromDataFrame(y9)
  y15gr<-makeGRangesFromDataFrame(y15)
  y22gr<-makeGRangesFromDataFrame(y22)
  
  alldmrs<-reduce(c(y9gr,y15gr,y22gr))
  mcols(alldmrs)$y9 = overlapsAny(alldmrs,y9gr)
  mcols(alldmrs)$y15 = overlapsAny(alldmrs,y15gr)
  mcols(alldmrs)$y22 = overlapsAny(alldmrs,y22gr)
  return(alldmrs)}

#alldmrs<-persistent_dmrs_fun(y9=y9_any2,y15=y15_any2,y22=y22_any2)
alldmrs<-persistent_dmrs_fun(y9=y9_dmrcate_combp,y15=y15_dmrcate_combp,y22=y22_dmrcate_combp)
persistent_dmrs<-alldmrs[alldmrs$y9==TRUE & alldmrs$y15 ==TRUE & alldmrs$y22==TRUE,]

#write.table(as.data.frame(alldmrs),file=paste0(OUTDIR,"/dmrs_from_combp_dmrcate_for_venn.txt"),sep="\t",row.names=F,quote=F)
write.table(as.data.frame(alldmrs),file=paste0(OUTDIR,"/dmrs_from_2methods_for_venn.txt"),sep="\t",row.names=F,quote=F)
#write.table(as.data.frame(persistent_dmrs),file=paste0(OUTDIR,"/persistent_dmrs_from_2methods_for_venn.txt"),sep="\t",row.names=F,quote=F)
write.table(as.data.frame(persistent_dmrs),file=paste0(OUTDIR,"/persistent_dmrs_from_combp_dmrcate_for_venn.txt"),sep="\t",row.names=F,quote=F)

any2_gr_fun <- function(any2) {
  any2_parsed<-any2
  any2_parsed$estimate<-ifelse(!is.na(any2_parsed$dmrcate_estimate),any2_parsed$dmrcate_estimate,ifelse(!is.na(any2_parsed$dmrff_estimate),any2_parsed$dmrff_estimate,any2_parsed$combp_estimate) )
  any2_parsed<-any2_parsed[,c("seqnames","start","end","estimate")]
  any2_parsed$estimate<-as.numeric(any2_parsed$estimate)
  any2_parsed_gr<-makeGRangesFromDataFrame(any2_parsed,keep.extra.columns = T)
}

y9_any2_gr<-any2_gr_fun(y9_dmrcate_combp)
y15_any2_gr<-any2_gr_fun(y15_dmrcate_combp)
y22_any2_gr<-any2_gr_fun(y22_dmrcate_combp)

any2_list<-list(y9=y9_any2_gr,y15=y15_any2_gr,y22=y22_any2_gr)
any2_combined<-make_dmr_master(any2_list)
any2_combined_df<-as.data.frame(any2_combined)
common_dmrs<-any2_combined_df[any2_combined_df$n_comparisons==3,]
write.table(as.data.frame(common_dmrs),file=paste0(OUTDIR,"/persistentdmrs_from_dmrcate_compbp_with_coordinates.txt"),sep="\t",row.names=F,quote=F)
write.table(as.data.frame(any2_combined_df),file=paste0(OUTDIR,"/dmrs_from_dmrcate_combp_with_coordinates.txt"),sep="\t",row.names=F,quote=F)

capture.output(
  print(sessionInfo()),
  file = paste0(OUTDIR, "/extract_dmrscalledbymultiplemethods_sessionInfo.txt"))

