#!/usr/bin/env Rscript

## Volcano plot of g:Profiler enrichment terms with REVIGO-based grouping
## Input:
##   - g:Profiler intersections CSV
##   - REVIGO tables for GO:BP, GO:CC, GO:MF
## Output:
##   - Volcano plot figure (PNG/PDF/SVG)
##   - Optional annotated tables for supplementary material

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggbreak)      # for scale_x_break
  library(ggnewscale)   # for new_scale_color
})

## -------------------------
## 1. Parameters & utilities
## -------------------------

# Significance threshold (adjustable)
sig_threshold_prof <- -log10(0.05)

# Function to compute odds ratio from contingency table
# a: intersection_size (in query & in term)
# term_size: total genes annotated with the term
# query_size: total genes in the query
# effective_domain_size: total genes in the background
calculate_OR <- function(a, term_size, query_size, effective_domain_size) {
  b <- term_size - a
  c <- query_size - a
  d <- effective_domain_size - term_size - query_size + a
  
  # Continuity correction
  if (b == 0 || c == 0 || d == 0) {
    b <- ifelse(b == 0, 0.5, b)
    c <- ifelse(c == 0, 0.5, c)
    d <- ifelse(d == 0, 0.5, d)
  }
  
  OR <- (a * d) / (b * c)
  return(OR)
}

# Color mapping by source (g:Profiler categories)
source_colors <- c(
  "GO:MF" = "#1b9e77",
  "GO:BP" = "#d95f02",
  "GO:CC" = "#7570b3",
  "KEGG"  = "#e7298a",
  "REAC"  = "#66a61e",
  "WP"    = "#e6ab02",
  "TF"    = "#a6761d",
  "MIRNA" = "#666666",
  "HPA"   = "#1f78b4",
  "CORUM" = "#b99f8a",
  "HP"    = "#fb9a99"
)
source_order <- names(source_colors)

# Paths
gprof_file <- "data/profiler/gProfiler_hsapiens_intersections.csv"

revigo_bp_file <- "data/profiler/GO_revigo_BP/Revigo_BP_Table.tsv"
revigo_cc_file <- "data/profiler/GO_revigo_CC/Revigo_CC_Table.tsv"
revigo_mf_file <- "data/profiler/GO_revigo_MF/Revigo_MF_Table.tsv"

output_dir <- "figures/tissue_enrichment"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## --------------------------------------
## 2. Build REVIGO group lookup (GO only)
## --------------------------------------

# Helper: process a single REVIGO table
process_revigo_table <- function(file_path) {
  read_tsv(file_path, show_col_types = FALSE) %>%
    mutate(
      Representative = na_if(Representative, "null"),
      Rep_full = if_else(
        !is.na(Representative),
        paste0("GO:", sprintf("%07d", as.numeric(Representative))),
        TermID
      )
    ) %>%
    # Map each TermID to its representative ID and representative name
    transmute(term_id = TermID, group_id = Rep_full) %>%
    left_join(
      read_tsv(file_path, show_col_types = FALSE) %>%
        select(TermID, Name),
      by = c("group_id" = "TermID")
    ) %>%
    transmute(term_id, group_id, group_id_name = Name)
}

revigo_bp_match <- process_revigo_table(revigo_bp_file)
revigo_cc_match <- process_revigo_table(revigo_cc_file)
revigo_mf_match <- process_revigo_table(revigo_mf_file)

revigo_key <- bind_rows(revigo_bp_match, revigo_cc_match, revigo_mf_match)

## ----------------------------------
## 3. Read g:Profiler and compute OR
## ----------------------------------

gProfiler_orig <- read_csv(gprof_file, show_col_types = FALSE)

# Ensure negative log10 adjusted p-values exist
if (!"negative_log10_of_adjusted_p_value" %in% colnames(gProfiler_orig)) {
  if ("adjusted_p_value" %in% colnames(gProfiler_orig)) {
    gProfiler_orig <- gProfiler_orig %>%
      mutate(
        negative_log10_of_adjusted_p_value = -log10(adjusted_p_value)
      )
  } else {
    stop("g:Profiler file must contain 'adjusted_p_value' or 'negative_log10_of_adjusted_p_value'.")
  }
}

gProfiler_OR <- gProfiler_orig %>%
  mutate(
    is_significant = negative_log10_of_adjusted_p_value > sig_threshold_prof,
    a = intersection_size,
    term_size = term_size,
    query_size = query_size,
    effective_domain_size = effective_domain_size,
    OR = calculate_OR(a, term_size, query_size, effective_domain_size)
  ) %>%
  # Attach REVIGO grouping (for GO terms)
  left_join(revigo_key, by = c("term_id")) %>%
  mutate(
    group_id = if_else(is.na(group_id), term_id, group_id),
    group_id_name = if_else(is.na(group_id_name), term_name, group_id_name),
    source = factor(source, levels = source_order),
    source_num = as.numeric(source),
    new_label = paste0("(", source_num, ") ", group_id_name),
    new_source_label = paste0("(", source_num, ") ", source)
  )

# Optional: save annotated table for Supplementary material
gProfiler_suppl_table <- gProfiler_orig %>%
  left_join(
    gProfiler_OR %>%
      select(source, term_id, group_id_name, OR),
    by = c("source", "term_id")
  )
# write.csv(gProfiler_suppl_table, file.path(output_dir, "gProfiler_enrichment_annotated.csv"), row.names = FALSE)

## --------------------------------------------
## 4. Select top terms per source for annotation
## --------------------------------------------

top_terms <- gProfiler_OR %>%
  filter(is_significant, !is.na(source)) %>%
  # one row per (source, OR, -log10p)
  group_by(source, OR, negative_log10_of_adjusted_p_value) %>%
  slice(1) %>%
  ungroup() %>%
  # one row per group_id_name
  arrange(desc(negative_log10_of_adjusted_p_value)) %>%
  group_by(group_id_name) %>%
  slice(1) %>%
  ungroup() %>%
  # top 5 per source
  arrange(desc(negative_log10_of_adjusted_p_value)) %>%
  group_by(source) %>%
  slice(1:5) %>%
  ungroup()

# Mark which terms are used in the plot annotation
gProfiler_OR_save <- gProfiler_OR %>%
  mutate(
    selected_for_graph = if_else(group_id_name %in% top_terms$group_id_name, 1L, 0L)
  )
# write.csv(gProfiler_OR_save, file.path(output_dir, "gProfiler_OR_all_terms.csv"), row.names = FALSE)

## -----------------------------------
## 5. Volcano plot (OR vs -log10 p)
## -----------------------------------

# Define label positions
y_start <- max(top_terms$negative_log10_of_adjusted_p_value, na.rm = TRUE) + 0.5
vertical_spacing <- 0.4

top_terms <- top_terms %>%
  arrange(desc(negative_log10_of_adjusted_p_value)) %>%
  mutate(
    label_x = 100,
    label_y = seq(from = y_start, by = -vertical_spacing, length.out = n())
  )

volcano_plot <- ggplot(gProfiler_OR,
                       aes(x = OR, y = negative_log10_of_adjusted_p_value)) +
  # Points: significant vs non-significant
  geom_point(aes(color = is_significant), size = 1) +
  geom_hline(yintercept = sig_threshold_prof, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "darkblue")) +
  theme(legend.position = "none") +
  # Optional x-axis break to de-emphasize extreme OR
  scale_x_break(c(300, 1200), ticklabels = c(1200, 1250), expand = TRUE) +
  xlim(0, 1250) +
  theme_classic() +
  labs(
    title = "Enrichment of functional / tissue terms (g:Profiler)",
    x = "Odds ratio (OR)",
    y = "-log10(adjusted p-value)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  # New color scale for annotations by source
  ggnewscale::new_scale_color() +
  geom_segment(
    data = top_terms,
    aes(
      x    = OR,
      xend = label_x,
      y    = negative_log10_of_adjusted_p_value,
      yend = label_y,
      color = source
    ),
    linetype = "dashed"
  ) +
  geom_text(
    data = top_terms,
    aes(x = label_x, y = label_y, label = group_id_name, color = source),
    size = 2.5,
    hjust = 0
  ) +
  scale_color_manual(values = source_colors)

print(volcano_plot)

ggsave(file.path(output_dir, "volcano_tissue_enrichment.png"),
       volcano_plot, width = 6, height = 4, dpi = 450, bg = "transparent")
ggsave(file.path(output_dir, "volcano_tissue_enrichment.pdf"),
       volcano_plot, width = 6, height = 6, device = "pdf", bg = "transparent")
ggsave(file.path(output_dir, "volcano_tissue_enrichment.svg"),
       volcano_plot, width = 6, height = 6, device = "svg", bg = "transparent")


