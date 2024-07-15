library(microbiomeMarker)

mm_lefse <- run_lefse(
  ps.prop.genus,
  wilcoxon_cutoff = 0.01,
  group = "BodySite",
  kw_cutoff = 0.01,
  multigrp_strat = TRUE,
  lda_cutoff = 4
)

mm_lefse


p_abd <- plot_abundance(mm_lefse, group = "BodySite")
p_abd

