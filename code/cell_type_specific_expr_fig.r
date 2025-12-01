library(data.table)
library(dplyr)
library(ggplot2)
library(scales)

# ---- Cell type order for plotting ----
cell_order <- c(
  "Pseudo-Bulk",
  "Neuronal",
  "Mast Cell",
  "Myeloid",
  "Adipocyte",
  "Lymphocyte",
  "VSMC",
  "Mural",
  "Pericyte",
  "Epicardial",
  "Endocardial",
  "Endothelial",
  "Lymphatic Endothelial",
  "Cardiac Endothelial",
  "Fibroblast",
  "Cardiomyocyte"
)

# ---- Load single-nucleus summary data ----
expr_file <- "/path/to/Reichart_LV_healthydonors__ExpressionValues_for_CellType_and_CellState__DruggabilityAllGenes.tsv"
# preprocessed from Pathogenic variants damage cell composition and single cell transcription in cardiomyopathies


df_whole <- fread(expr_file) %>%
  filter(!is.na(gene_ids)) %>%
  # keep only prioritized genes
  filter(gene_ids %in% summary_table_fin_drugnome_ai$genes_highest_score)

# Rename convenience columns
df_whole <- df_whole %>%
  mutate(
    Gene       = gene_ids,
    Percentage = Percent_nonZero_Expression,
    Expression = Mean_logNormalized_Expression
  )

# Optionally rescale expression per gene over [0, 1]
for (gene in unique(df_whole$Gene)) {
  idx <- df_whole$Gene == gene
  if (!all(is.na(df_whole$Expression[idx]))) {
    df_whole$Expression[idx] <- scales::rescale(
      df_whole$Expression[idx],
      to = c(0, 1)
    )
  }
}

# ---- Prepare data frame for plotting ----
df <- df_whole %>%
  group_by(gene_ids) %>%
  mutate(Expression_scaled = scales::rescale(Mean_logNormalized_Expression, to = c(0, 1))) %>%
  ungroup()

# Order cell types
df$CellType <- factor(df$CellType, levels = cell_order)

# Order genes by mean percent expression across all cell types
gene_order <- df %>%
  group_by(gene_ids) %>%
  summarize(mean_pct = mean(Percentage, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_pct)) %>%
  pull(gene_ids)

# y-axis: genes from top (highest mean_pct) to bottom
df$gene_ids <- factor(df$gene_ids, levels = rev(gene_order))

# ---- Dot plot ----
p <- ggplot(df, aes(x = CellType, y = gene_ids)) +
  # light grid tiles
  geom_tile(
    fill  = "white",
    color = "grey90",
    width = 0.95,
    height = 0.95
  ) +
  # dot: size = % expressing; fill = scaled expression
  geom_point(
    aes(size = Percentage, fill = Expression_scaled),
    shape = 21,
    stroke = 0.4,
    color  = "black"
  ) +
  scale_fill_gradient2(
    low  = "blue",
    mid  = "white",
    high = "red",
    midpoint = 0.5,
    limits   = c(0, 1),
    name     = "Scaled\nmean\nexpression"
  ) +
  scale_size_continuous(
    range  = c(0, 7),
    limits = c(0, 100),
    name   = "% nuclei\nexpressing"
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      vjust = 1,
      size  = 11
    ),
    axis.text.y = element_text(size = 11),
    legend.position   = "right",
    panel.grid        = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA)
  ) +
  labs(
    x = "Cell type",
    y = "Gene",
    title = ""
  )

# Print or save
print(p)
# ggsave("ExtendedData_Fig8_celltype_dots.pdf", p, width = 8, height = 10)