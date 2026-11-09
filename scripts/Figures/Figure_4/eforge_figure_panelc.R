
library(ggplot2)
library(dplyr)
library(Cairo)

setwd("path to working directory")

# ------------------------------------------------------------
# 1. Read eFORGE results
# ------------------------------------------------------------

all_results <- read.table(
  "Unnamed.850k.erc2-chromatin15state-all.chart.tsv.gz",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Filter and prepare data
# ------------------------------------------------------------

tissue_order <- c(
  "Blood", "Heart", "ES-derived", "ESC",
  "iPSC", "Mesenchymal", "Spleen", "Thymus"
)

sig_results <- all_results %>%
  filter(
    Qvalue < 0.01,
    Tissue %in% tissue_order,
    !is.na(Pvalue),
    Pvalue > 0
  ) %>%
  mutate(
    logP = -log10(Pvalue),
    Tissue = factor(Tissue, levels = tissue_order)
  ) %>%
  arrange(Tissue, desc(logP)) %>%
  mutate(
    Tissue_Cell_num = row_number(),
    CellLabel = gsub(
      "hematopoietic stem cells",
      "HSC",
      Cell
    )
  )


# ------------------------------------------------------------
# 3. Create background blocks for tissue groups
# ------------------------------------------------------------

bg_blocks <- sig_results %>%
  group_by(Tissue) %>%
  summarise(
    xmin = min(Tissue_Cell_num) - 0.5,
    xmax = max(Tissue_Cell_num) + 0.5,
    .groups = "drop"
  )

# Position tissue labels within their background blocks

x_min <- min(sig_results$logP, na.rm = TRUE)
x_max <- max(sig_results$logP, na.rm = TRUE)

bg_blocks <- bg_blocks %>%
  mutate(
    label_y = x_min + 0.40 * (x_max - x_min),
    label_x = xmax - 0.8
  )

# ------------------------------------------------------------
# 4. Chromatin-state legend
# ------------------------------------------------------------

state_order <- c("Enh", "EnhBiv", "TssAFlnk")

sig_results <- sig_results %>%
  mutate(
    Datatype = factor(Datatype, levels = state_order)
  )

# ------------------------------------------------------------
# 5. Plot
# ------------------------------------------------------------

p <- ggplot(
  sig_results,
  aes(x = Tissue_Cell_num, y = logP)
) +
  
  # Tissue background bands
  geom_rect(
    data = bg_blocks,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = -Inf,
      ymax = Inf,
      fill = Tissue
    ),
    alpha = 0.15,
    inherit.aes = FALSE
  ) +
  
  # Enrichment points
  geom_point(
    aes(
      color = Tissue,
      shape = Datatype
    ),
    size = 1.6
  ) +
  
  
  scale_shape_discrete(
    breaks = state_order,
    labels = state_order
  ) +
  # Cell descriptions
  scale_x_continuous(
    breaks = sig_results$Tissue_Cell_num,
    labels = sig_results$CellLabel,
    expand = expansion(add = 0.5)
  ) +
  
  coord_flip() +
  
  labs(
    x = NULL,
    y = expression(-"log"[10]~P),
    shape = NULL
  ) +
  
  guides(
    color = "none",
    fill = "none",
    shape = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 2.2,
        color = "black"
      )
    )
  ) +
  
  theme_classic(
    base_family = "Arial",
    base_size = 8
  ) +
  
  theme(
    # Cell descriptions
    axis.text.y = element_text(
      size = 6,
      family = "Arial",
      margin = margin(r = 4, unit = "pt")
    ),
    
    # Horizontal significance-axis numbers
    axis.text.x = element_text(
      size = 7,
      family = "Arial"
    ),
    
    # Significance-axis title
    axis.title.x = element_text(
      size = 8,
      family = "Arial",
      margin = margin(t = 4, unit = "pt")
    ),
    
    axis.title.y = element_blank(),
    
    # Chromatin-state legend
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 8,
      family = "Arial",
      face = "plain"
    ),
    
    legend.key.width = grid::unit(2, "mm"),
    legend.key.height = grid::unit(4, "mm"),
    legend.spacing.x = grid::unit(8, "mm"),
    
    legend.margin = margin(
      5, 4, 4, 4,
      unit = "pt"
    ),
    
    plot.margin = margin(
      5, 5, 5, 5,
      unit = "pt"
    )
  )

# ------------------------------------------------------------
# 6. Add tissue group labels
# ------------------------------------------------------------

p <- p +
  geom_text(
    data = bg_blocks,
    aes(
      x = label_x,
      y = label_y,
      label = Tissue
    ),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 0.5,
    size = 6,
    size.unit = "pt",
    family = "Arial",
    fontface = "bold",
    color = "black"
  )

# ------------------------------------------------------------
# 7. Export PDF
# ------------------------------------------------------------

CairoPDF(
  "eFORGE_y22_sig_results_panelc.pdf",
  family = "Arial",
  width = 5.5,
  height = 4.75
)

library(patchwork)


tag_theme<-theme(plot.tag=element_text(
  face="bold",family="Arial",size=14))
# Assemble plots
pc <- (p+tag_theme)

p+ plot_annotation(tag_levels=list('c'))



print(pc)
dev.off()
