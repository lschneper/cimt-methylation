library(readxl)
library(ggplot2)
library(forcats)
library(dplyr)
library(ComplexHeatmap)
library(UpSetR)
library(tidyr)
library(simplifyEnrichment)
library(grid)
library(ggtext)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)
library(tidyr)
library(Cairo)

# DMRs 
INDIR="/Users/monet/Desktop/Eve_2026/Sept2026/Fig4/fig4/"
OUTDIR="/Users/monet/Desktop/Eve_2026/Sept2026/Fig6_test"
dir.create(OUTDIR)
####################################################################
##### Read in data ##################################
####################################################################
setwd(INDIR)
y9DMR<-readRDS("y9_DMRs_DMRcate_combp_chipseeker.RDS")
y15DMR<-readRDS("y15_DMRs_DMRcate_combp_chipseeker.RDS")
y22DMR<-readRDS("y22_DMRs_DMRcate_combp_chipseeker.RDS")

y9EWAS<-readRDS("y9_meta_sensitivity_model1_withRE_metalistdf_output_bacon_chipseeker.RDS")
y15EWAS<-readRDS("y15_meta_sensitivity_model1_withRE_metalistdf_output_bacon_chipseeker.RDS")
y22EWAS<-readRDS("y22_epic_age_sex_7svs_newmiddleton_031426.ewasoutput_man_annotated_chipseeker_chipseeker.RDS")
#######################################################################
#### Functions #######################################################
#######################################################################
# Convert to ENTREZID
convert_genes <- function(x) {
  bitr(
    x,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID) %>%
    pull(ENTREZID)
}

# get DMR promoters & genes
get_promoter_gene_genes <- function(x) {
  
  x %>%
    dplyr::select(
      overlapping_promoterIDs,
      overlapping_geneIDs
    ) %>%
    tidyr::pivot_longer(
      cols = everything(),
      values_to = "gene"
    ) %>%
    dplyr::filter(!is.na(gene)) %>%
    tidyr::separate_rows(gene, sep = ",") %>%
    dplyr::mutate(
      gene = trimws(gene),
      gene = sub("\\.[0-9]+$", "", gene)
    ) %>%
    dplyr::filter(gene != "") %>%
    dplyr::distinct(gene) %>%
    dplyr::pull(gene)
}
# GO Biological Process and Molecular Function
run_go <- function(genes, universe, ontology = "BP") {
  
  clusterProfiler::enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
}

run_reactome <- function(genes, universe) {
  ReactomePA::enrichPathway(
    gene          = genes,
    universe      = universe,
    organism      = "human",
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    qvalueCutoff  = 1,
    readable      = TRUE
  )
}


#################################################################
# Extract overlapping promoters & genes in DMRs and EWAS universe
#################################################################

dmr_promoter_gene <- list(
  y9  = get_promoter_gene_genes(y9DMR),
  y15 = get_promoter_gene_genes(y15DMR),
  y22 = get_promoter_gene_genes(y22DMR)
)

dmr_promoter_gene_entrez <- lapply(
  dmr_promoter_gene,
  convert_genes
)
ewas_promoter_gene_universe <- list(
  y9  = get_promoter_gene_genes(y9EWAS),
  y15 = get_promoter_gene_genes(y15EWAS),
  y22 = get_promoter_gene_genes(y22EWAS)
)
ewas_promoter_gene_universe_ENTREZ <-
  lapply(ewas_promoter_gene_universe, convert_genes)
mapply(
  function(fg, bg) length(setdiff(fg, bg)),
  dmr_promoter_gene_entrez,
  ewas_promoter_gene_universe_ENTREZ
)

######################################################################
#### BP #################
######################################################################
BP_gp <- list(
  y9 = run_go(
    dmr_promoter_gene_entrez$y9,
    ewas_promoter_gene_universe_ENTREZ$y9,
    "BP"
  ),
  
  y15 = run_go(
    dmr_promoter_gene_entrez$y15,
    ewas_promoter_gene_universe_ENTREZ$y15,
    "BP"
  ),
  
  y22 = run_go(
    dmr_promoter_gene_entrez$y22,
    ewas_promoter_gene_universe_ENTREZ$y22,
    "BP"
  )
)

BP_gp_all <- dplyr::bind_rows(
  lapply(names(BP_gp), function(x) {
    as.data.frame(BP_gp[[x]]) %>%
      dplyr::mutate(Age = x)
  })
)

BP_for_simplify <- BP_gp_all %>%
  dplyr::filter(
    pvalue < 0.05,
    Count >= 2
  ) %>%
  split(.$Age) %>%
  lapply(function(x) unique(x$ID))

sapply(BP_for_simplify, length)
write.xlsx(BP_gp_all,file=paste0(OUTDIR,"/BP.xlsx"))
BP_clusters <- simplifyEnrichment::simplifyGOFromMultipleLists(
  BP_for_simplify,
  ont = "BP",
  db = "org.Hs.eg.db",
  method = "binary_cut",
  column_title = "GO Biological Process"
)

BP_cluster_info <- BP_clusters %>%
  dplyr::rename(ID = id) %>%
  dplyr::left_join(
    BP_gp_all %>%
      dplyr::filter(
        pvalue < 0.05,
        Count >= 2
      ) %>%
      dplyr::group_by(ID, Description) %>%
      dplyr::summarise(
        n_ages = dplyr::n_distinct(Age),
        ages = paste(sort(unique(Age)), collapse = ", "),
        min_p = min(pvalue),
        min_FDR = min(p.adjust),
        max_count = max(Count),
        .groups = "drop"
      ),
    by = "ID"
  ) %>%
  dplyr::arrange(
    cluster,
    dplyr::desc(n_ages),
    min_p
  )
as.data.frame(BP_cluster_info)

BP_representatives <- BP_cluster_info %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()
BP_representatives
as.data.frame(BP_representatives)

selected_BP <- BP_representatives$ID
## Preserve the cluster order from simplifyEnrichment
BP_order <- BP_representatives %>%
  arrange(cluster) %>%
  pull(Description)

## Pull results for those terms from the FULL enrichment results
BP_plot <- BP_gp_all %>%
  filter(ID %in% selected_BP) %>%
  mutate(
    Age = factor(
      Age,
      levels = c("y9", "y15", "y22"),
      labels = c("Age 9", "Age 15", "Age 22")
    ),
    score = -log10(pvalue),
    FDR_sig = p.adjust < 0.05
  )

BP_plot %>%
  count(Description) %>%
  arrange(n)

expand_grid(
  ID = selected_BP,
  Age = factor(
    c("Age 9", "Age 15", "Age 22"),
    levels = c("Age 9", "Age 15", "Age 22")
  )
) %>%
  anti_join(
    BP_plot %>% dplyr::select(ID, Age),
    by = c("ID", "Age")
  )
missing_BP <- expand_grid(
  ID = selected_BP,
  Age = factor(
    c("Age 9", "Age 15", "Age 22"),
    levels = c("Age 9", "Age 15", "Age 22")
  )
) %>%
  anti_join(
    BP_plot %>% dplyr::select(ID, Age),
    by = c("ID", "Age")
  )

missing_BP
#### recode the missing values
BP_labels <- BP_representatives %>%
  dplyr::select(ID, Description, cluster)

BP_plot_complete <- tidyr::expand_grid(
  ID = selected_BP,
  Age = factor(
    c("Age 9", "Age 15", "Age 22"),
    levels = c("Age 9", "Age 15", "Age 22")
  )
) %>%
  dplyr::left_join(
    BP_plot %>%
      dplyr::select(
        ID, Age, pvalue, p.adjust,
        Count, score, FDR_sig
      ),
    by = c("ID", "Age")
  ) %>%
  dplyr::left_join(
    BP_labels,
    by = "ID"
  )

table(is.na(BP_plot_complete$pvalue))
BP_plot_complete <- BP_plot_complete %>%
  dplyr::mutate(
    Description = factor(
      Description,
      levels = rev(BP_representatives$Description)
    )
  )
BP_plot_complete <- BP_plot_complete %>%
  dplyr::mutate(
    cell_label = dplyr::case_when(
      is.na(Count) ~ "",
      !is.na(FDR_sig) & FDR_sig ~ paste0(Count, "*"),
      TRUE ~ as.character(Count)
    ),
    label_color = dplyr::if_else(
      !is.na(score) & score >= 3,
      "white",
      "black"
    )
  )



#######################################################################
#### MF #################
#######################################################################
MF_gp <- list(
  y9 = run_go(
    dmr_promoter_gene_entrez$y9,
    ewas_promoter_gene_universe_ENTREZ$y9,
    "MF"
  ),
  
  y15 = run_go(
    dmr_promoter_gene_entrez$y15,
    ewas_promoter_gene_universe_ENTREZ$y15,
    "MF"
  ),
  
  y22 = run_go(
    dmr_promoter_gene_entrez$y22,
    ewas_promoter_gene_universe_ENTREZ$y22,
    "MF"
  )
)
MF_gp_all <- dplyr::bind_rows(
  lapply(names(MF_gp), function(x) {
    as.data.frame(MF_gp[[x]]) %>%
      dplyr::mutate(Age = x)
  })
)

MF_for_simplify <- MF_gp_all %>%
  dplyr::filter(
    pvalue < 0.05,
    Count >= 2
  ) %>%
  split(.$Age) %>%
  lapply(function(x) unique(x$ID))

sapply(MF_for_simplify, length)
write.xlsx(MF_gp,file=paste0(OUTDIR,"/MF.xlsx"))

MF_ids <- unique(unlist(MF_for_simplify))

length(MF_ids)
# should be 35
set.seed(123)
MF_sim <- simplifyEnrichment::GO_similarity(
  MF_ids,
  ont = "MF",
  db = "org.Hs.eg.db"
)

set.seed(123)

MF_clusters2 <- simplifyEnrichment::simplifyGO(
  MF_sim,
  method = "binary_cut",
  verbose = TRUE
)


########
# Make the figure

MF_cluster_info <- MF_clusters2 %>%
  dplyr::rename(ID = id) %>%
  dplyr::left_join(
    MF_gp_all %>%
      dplyr::filter(
        pvalue < 0.05,
        Count >= 2
      ) %>%
      dplyr::group_by(ID, Description) %>%
      dplyr::summarise(
        n_ages = dplyr::n_distinct(Age),
        ages = paste(sort(unique(Age)), collapse = ", "),
        min_p = min(pvalue),
        min_FDR = min(p.adjust),
        max_count = max(Count),
        .groups = "drop"
      ),
    by = "ID"
  ) %>%
  dplyr::arrange(
    cluster,
    dplyr::desc(n_ages),
    min_p
  )


as.data.frame(MF_cluster_info)
table(MF_cluster_info$cluster)
MF_representatives <- MF_cluster_info %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()
MF_representatives
as.data.frame(MF_representatives)

## IDs selected as representatives of the semantic clusters
selected_MF <- MF_representatives$ID

## Preserve the cluster order from simplifyEnrichment
MF_order <- MF_representatives %>%
  arrange(cluster) %>%
  pull(Description)

## Pull results for those terms from the FULL enrichment results
MF_plot <- MF_gp_all %>%
  filter(ID %in% selected_MF) %>%
  mutate(
    Age = factor(
      Age,
      levels = c("y9", "y15", "y22"),
      labels = c("Age 9", "Age 15", "Age 22")
    ),
    score = -log10(pvalue),
    FDR_sig = p.adjust < 0.05
  )

MF_plot %>%
  count(Description) %>%
  arrange(n)

expand_grid(
  ID = selected_MF,
  Age = factor(
    c("Age 9", "Age 15", "Age 22"),
    levels = c("Age 9", "Age 15", "Age 22")
  )
) %>%
  anti_join(
    MF_plot %>% dplyr::select(ID, Age),
    by = c("ID", "Age")
  )
missing_MF <- expand_grid(
  ID = selected_MF,
  Age = factor(
    c("Age 9", "Age 15", "Age 22"),
    levels = c("Age 9", "Age 15", "Age 22")
  )
) %>%
  anti_join(
    MF_plot %>% dplyr::select(ID, Age),
    by = c("ID", "Age")
  )

missing_MF
#### recode the missing values
MF_labels <- MF_representatives %>%
  dplyr::select(ID, Description, cluster)

MF_plot_complete <- tidyr::expand_grid(
  ID = selected_MF,
  Age = factor(
    c("Age 9", "Age 15", "Age 22"),
    levels = c("Age 9", "Age 15", "Age 22")
  )
) %>%
  dplyr::left_join(
    MF_plot %>%
      dplyr::select(
        ID, Age, pvalue, p.adjust,
        Count, score, FDR_sig
      ),
    by = c("ID", "Age")
  ) %>%
  dplyr::left_join(
    MF_labels,
    by = "ID"
  )

table(is.na(MF_plot_complete$pvalue))
# Plot
MF_plot_complete <- MF_plot_complete %>%
  dplyr::mutate(
    Description = factor(
      Description,
      levels = rev(MF_representatives$Description)
    )
  )

# Set text color to white if score > 3
MF_plot_complete <- MF_plot_complete %>%
  dplyr::mutate(
    cell_label = dplyr::case_when(
      is.na(Count) ~ "",
      !is.na(FDR_sig) & FDR_sig ~ paste0(Count, "*"),
      TRUE ~ as.character(Count)
    ),
    label_color = dplyr::if_else(
      !is.na(score) & score >= 3,
      "white",
      "black"
    )
  )

MF_cluster_info %>%
  dplyr::filter(
    grepl(
      "transporter regulator|ion channel regulator",
      Description,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::select(
    cluster, ID, Description,
    n_ages, ages, min_p, min_FDR
  )


########################################################################
### REACTOME ##############
########################################################################

Reactome_gp <- list(
  y9 = run_reactome(
    dmr_promoter_gene_entrez$y9,
    ewas_promoter_gene_universe_ENTREZ$y9
  ),
  y15 = run_reactome(
    dmr_promoter_gene_entrez$y15,
    ewas_promoter_gene_universe_ENTREZ$y15
  ),
  y22 = run_reactome(
    dmr_promoter_gene_entrez$y22,
    ewas_promoter_gene_universe_ENTREZ$y22
  )
)

Reactome_gp_all <- dplyr::bind_rows(
  lapply(names(Reactome_gp), function(x) {
    as.data.frame(Reactome_gp[[x]]) %>%
      dplyr::mutate(Age = x)
  })
)

write.xlsx(Reactome_gp_all,file=paste0(OUTDIR,"/Reactome.xlsx"))

Reactome_gp_all %>%
  dplyr::group_by(Age) %>%
  dplyr::summarise(
    n_tested = dplyr::n(),
    min_p = min(pvalue),
    min_FDR = min(p.adjust),
    n_p05 = sum(pvalue < 0.05),
    n_FDR05 = sum(p.adjust < 0.05),
    .groups = "drop"
  )

Reactome_gp_candidates <- Reactome_gp_all %>%
  dplyr::filter(
    pvalue < 0.05,
    Count >= 2
  )

Reactome_gp_candidates %>%
  dplyr::arrange(Age, pvalue) %>%
  dplyr::select(
    Age, Description, Count,
    pvalue, p.adjust, geneID
  )

Reactome_gp_term_summary <- Reactome_gp_candidates %>%
  dplyr::group_by(ID, Description) %>%
  dplyr::summarise(
    n_ages = dplyr::n_distinct(Age),
    ages = paste(sort(unique(Age)), collapse = ", "),
    min_p = min(pvalue),
    min_FDR = min(p.adjust),
    max_count = max(Count),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_ages),
    min_p
  )

react_candidates <- Reactome_gp_all %>%
  dplyr::filter(
    pvalue < 0.05,
    Count >= 2
  )
react_labels <- react_candidates %>%
  dplyr::distinct(ID, Description)

library(reactome.db)
react_ids <- unique(react_candidates$ID)

react_map <- AnnotationDbi::select(
  reactome.db,
  keys = react_ids,
  keytype = "PATHID",
  columns = c(
    "PATHID",
    "PATHNAME",
    "ENTREZID"
  )
) %>%
  dplyr::filter(!is.na(ENTREZID))

# Check coverage
setdiff(
  react_ids,
  unique(react_map$PATHID)
)

# Cereate Reactome gene sets
react_sets <- split(
  react_map$ENTREZID,
  react_map$PATHID
)

react_sets <- lapply(react_sets, unique)

length(react_sets)
sapply(react_sets, length) |> summary()


# Use Jaccard similarity to calculate gene member similarity matrix
jaccard <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}

react_names <- names(react_sets)

react_sim <- outer(
  react_names,
  react_names,
  Vectorize(function(x, y) {
    jaccard(react_sets[[x]], react_sets[[y]])
  })
)

dimnames(react_sim) <- list(react_names, react_names)

#### Make table of unique pairs
react_pairs <- as.data.frame(as.table(react_sim)) %>%
  dplyr::rename(
    ID1 = Var1,
    ID2 = Var2,
    jaccard = Freq
  ) %>%
  dplyr::mutate(
    ID1 = as.character(ID1),
    ID2 = as.character(ID2)
  ) %>%
  dplyr::filter(ID1 < ID2) %>%
  dplyr::left_join(
    react_labels %>%
      dplyr::rename(
        ID1 = ID,
        Description1 = Description
      ),
    by = "ID1"
  ) %>%
  dplyr::left_join(
    react_labels %>%
      dplyr::rename(
        ID2 = ID,
        Description2 = Description
      ),
    by = "ID2"
  ) %>%
  dplyr::arrange(dplyr::desc(jaccard))

# Calculate containment
react_pairs <- react_pairs %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    n1 = length(react_sets[[ID1]]),
    n2 = length(react_sets[[ID2]]),
    n_overlap = length(
      intersect(
        react_sets[[ID1]],
        react_sets[[ID2]]
      )
    ),
    containment1 = n_overlap / n1,
    containment2 = n_overlap / n2,
    max_containment = max(containment1, containment2),
    min_containment = min(containment1, containment2)
  ) %>%
  dplyr::ungroup()

react_pairs %>%
  dplyr::arrange(dplyr::desc(jaccard)) %>%
  dplyr::select(
    Description1, Description2,
    n1, n2, n_overlap,
    jaccard, max_containment
  ) %>%
  head(30)

### Test collapsing
react_edges <- react_pairs %>%
  dplyr::filter(
    max_containment >= 0.8,
    jaccard >= 0.5
  ) %>%
  dplyr::select(ID1, ID2)
library(igraph)

g <- igraph::graph_from_data_frame(
  react_edges,
  directed = FALSE,
  vertices = react_ids
)

react_membership <- igraph::components(g)$membership

react_clusters <- tibble::tibble(
  ID = names(react_membership),
  cluster = as.integer(react_membership)
)

table(react_clusters$cluster)
# Join enrichment result
react_cluster_info <- react_clusters %>%
  dplyr::left_join(
    react_candidates %>%
      dplyr::group_by(ID, Description) %>%
      dplyr::summarise(
        n_ages = dplyr::n_distinct(Age),
        ages = paste(sort(unique(Age)), collapse = ", "),
        min_p = min(pvalue),
        min_FDR = min(p.adjust),
        max_count = max(Count),
        .groups = "drop"
      ),
    by = "ID"
  ) %>%
  dplyr::arrange(
    cluster,
    dplyr::desc(n_ages),
    min_p
  )
##### 
react_cluster_info %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(cluster_size = dplyr::n()) %>%
  dplyr::filter(cluster_size > 1) %>%
  dplyr::select(
    cluster, ID, Description,
    n_ages, ages, min_p, min_FDR
  ) %>%
  as.data.frame()

react_representatives <- react_cluster_info %>%
  dplyr::arrange(
    cluster,
    dplyr::desc(n_ages),
    min_p,
    Description
  ) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()
selected_react <- react_representatives$ID
react_order <- react_representatives %>%
  dplyr::arrange(cluster) %>%
  dplyr::pull(Description)


react_plot <- Reactome_gp_all %>%
  dplyr::filter(ID %in% selected_react) %>%
  dplyr::mutate(
    Age = factor(
      Age,
      levels = c("y9", "y15", "y22"),
      labels = c("Age 9", "Age 15", "Age 22")
    ),
    score = -log10(pvalue),
    FDR_sig = p.adjust < 0.05
  )
length(selected_react)
react_labels <- react_representatives %>%
  dplyr::select(ID, Description, cluster)

react_plot_complete <- tidyr::expand_grid(
  ID = selected_react,
  Age = factor(
    c("Age 9", "Age 15", "Age 22"),
    levels = c("Age 9", "Age 15", "Age 22")
  )
) %>%
  dplyr::left_join(
    react_plot %>%
      dplyr::select(
        ID, Age, pvalue, p.adjust,
        Count, score, FDR_sig
      ),
    by = c("ID", "Age")
  ) %>%
  dplyr::left_join(
    react_labels,
    by = "ID"
  )
nrow(react_plot_complete)
# 51
sum(is.na(react_plot_complete$score))
# 11

react_plot_complete <- react_plot_complete %>%
  dplyr::mutate(
    cell_label = dplyr::case_when(
      is.na(Count) ~ "",
      FDR_sig ~ paste0(Count, "*"),
      TRUE ~ as.character(Count)
    ),
    label_color = dplyr::if_else(
      !is.na(score) & score >= 3,
      "white",
      "black"
    ),
    label_font = dplyr::if_else(
      !is.na(FDR_sig) & FDR_sig,
      "bold",
      "plain"
    ),
    Description = factor(
      Description,
      levels = rev(react_representatives$Description)
    )
  )

react_plot_complete %>%
  dplyr::filter(is.na(score)) %>%
  dplyr::select(Description, Age)
missing_check <- react_plot_complete %>%
  dplyr::filter(is.na(score)) %>%
  dplyr::mutate(
    age_key = dplyr::recode(
      as.character(Age),
      "Age 9"  = "y9",
      "Age 15" = "y15",
      "Age 22" = "y22"
    )
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    n_fg_genes = length(
      intersect(
        react_sets[[ID]],
        dmr_promoter_gene_entrez[[age_key]]
      )
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(Description, Age, n_fg_genes)

missing_check




### IPA
files<-list.files(pattern="IPA")
files<-files[grep("pdf",files,invert=T)]
samples<-c("Age 15","Age 22","Age 9")
filelist<-list()
for(i in 1:length(files)){
  tmp<-as.data.frame(read_excel(files[i],skip=1))
  tmp$Count<-lengths(strsplit(tmp$Molecules,","))
  tmp$label_color<-dplyr::if_else(
    !is.na(tmp$`-log(p-value)`) & tmp$`-log(p-value)` >= 3,
    "white",
    "black"
  )
  
  filelist[[i]]<-tmp
}

names(filelist)<-samples
ipaproc<-vector()
for(i in 1:length(filelist)){
  ipaproc<-c(ipaproc,filelist[[i]][filelist[[i]]$`-log(p-value)`>=1.3, "Ingenuity Canonical Pathways"])}
ipaproc<-ipaproc[!duplicated(ipaproc)]
ipaproc<-gsub("ipa_","",ipaproc)
ipadf<-data.frame(ipa_Term=ipaproc)

filelistfiltered<-list()
for(i in 1:length(filelist)){
  filelistfiltered[[i]]<-filelist[[i]][filelist[[i]]$`Ingenuity Canonical Pathways` %in% ipaproc,]
  filelistfiltered[[i]]$Age<-rep(names(filelist)[i],nrow(filelistfiltered[[i]]))
}

filelistfiltereddf<-do.call(rbind.data.frame,filelistfiltered)
filelistfiltereddf$Age<-factor(filelistfiltereddf$Age,levels=c("Age 9","Age 15","Age 22"))
filelistfiltereddfo<-filelistfiltereddf[order(filelistfiltereddf$`Ingenuity Canonical Pathways`,filelistfiltereddf$Age),]
###########
top_IPA <- lapply(filelist, function(x) {
  x %>%
    filter(!is.na(`-log(p-value)`)) %>%
    slice_max(
      order_by = `-log(p-value)`,
      n = 7,
      with_ties = FALSE
    )
})

top_terms <- unique(
  unlist(
    lapply(top_IPA, function(x)
      x$`Ingenuity Canonical Pathways`
    )
  )
)

length(top_terms)
IPA_plot_list <- list()

for(i in seq_along(filelist)) {
  
  IPA_plot_list[[i]] <- filelist[[i]] %>%
    filter(`Ingenuity Canonical Pathways` %in% top_terms) %>%
    mutate(Age = names(filelist)[i])
}

IPA_plot <- bind_rows(IPA_plot_list)

IPA_plot$Age <- factor(
  IPA_plot$Age,
  levels = c("Age 9", "Age 15", "Age 22")
)

sapply(c(5, 7, 10, 15), function(n) {
  
  terms <- unique(
    unlist(
      lapply(filelist, function(x) {
        x %>%
          filter(!is.na(`-log(p-value)`)) %>%
          slice_max(
            order_by = `-log(p-value)`,
            n = n,
            with_ties = FALSE
          ) %>%
          pull(`Ingenuity Canonical Pathways`)
      })
    )
  )
  
  length(terms)
})





IPA_sig <- bind_rows(
  lapply(seq_along(filelist), function(i) {
    filelist[[i]] %>%
      mutate(Age = names(filelist)[i])
  })
) %>%
  filter(`-log(p-value)` >= 1.3)

IPA_summary <- IPA_sig %>%
  group_by(`Ingenuity Canonical Pathways`) %>%
  summarise(
    n_ages = n_distinct(Age),
    max_score = max(`-log(p-value)`, na.rm = TRUE),
    ages = paste(Age, collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_ages), desc(max_score))

keep_IPA <- IPA_summary %>%
  filter(n_ages >= 2) %>%
  pull(`Ingenuity Canonical Pathways`)

length(keep_IPA)
# 20
IPA_plot <- bind_rows(
  lapply(seq_along(filelist), function(i) {
    filelist[[i]] %>%
      filter(`Ingenuity Canonical Pathways` %in% keep_IPA) %>%
      mutate(Age = names(filelist)[i])
  })
)

IPA_plot$Age <- factor(
  IPA_plot$Age,
  levels = c("Age 9", "Age 15", "Age 22")
)


IPA_plot <- IPA_plot %>%
  tidyr::complete(
    `Ingenuity Canonical Pathways` = keep_IPA,
    Age = factor(
      c("Age 9", "Age 15", "Age 22"),
      levels = c("Age 9", "Age 15", "Age 22")
    )
  )



IPA_plot <- IPA_plot %>%
  mutate(
    cell_label = ifelse(
      is.na(`-log(p-value)`),
      "",
      as.character(Count)
    ),
    label_color = ifelse(
      !is.na(`-log(p-value)`) & `-log(p-value)` >= 3,
      "white",
      "black"
    )
  )

IPA_order <- IPA_summary %>%
  filter(n_ages >= 2) %>%
  arrange(desc(n_ages), desc(max_score)) %>%
  pull(`Ingenuity Canonical Pathways`)

IPA_plot$`Ingenuity Canonical Pathways` <- factor(
  IPA_plot$`Ingenuity Canonical Pathways`,
  levels = rev(IPA_order)
)



sourcedata <- list(
  BP = BP_plot_complete,
  MF = MF_plot_complete,
  Reactome = react_plot_complete,
  ipa = IPA_plot
)

openxlsx::write.xlsx(
  sourcedata,
  file = file.path(OUTDIR, "Fig6_sourcedata.xlsx"),
  overwrite = TRUE
)




writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "Fig6.sessionInfo.txt")
)
