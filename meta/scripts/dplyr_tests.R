# NOTE: The following files are large-ish and are so meta/data/ is in .gitignore
# IDH2mut_v_NBM_classification_comparison.bed
# 2607_mc_hmc_perc_meth_test.txt
# IDH2mut_v_NBM_test.txt

devtools::load_all()
library(ggplot2)

################################################################################
# Example of regions with continuous data and no strand (fairly small ~164K rows)
# Runs read to vis in ~15s
# Runs read to annotate in ~7s

file = '../data/IDH2mut_v_NBM_test.txt'
d = read_bed(file = file, genome = 'hg19', stranded = FALSE, use.score = TRUE)
annotations = c('cpgs', 'detailed_genes')
i = intersect_annotations(
  regions = d,
  annotations = annotations,
  genome = 'hg19',
  ignore.strand = T)
t = annotate_intersections(
  regions = d,
  intersections = i,
  use.score = T)

aggregation_over_features = dplyr::summarize(
  dplyr::group_by(t, annot_type, annot_id),
  mean = mean(data_score),
  sd = sd(data_score))

plot =
  ggplot(aggregation_over_features, aes(mean)) +
  geom_histogram(binwidth=5, aes(y=..density..)) +
  facet_wrap(~annot_type) +
  theme_bw()
ggsave(file = '../figures/IDH2mut_v_NBM_DMdiff_over_annots.png',
  plot = plot, width = 12, height = 12, dpi = 300)

################################################################################
# Example of regions with continuous data and strand (fairly large, ~4M rows)
# Runs read to vis in ~2m
# Runs read to annotate in ~1:53m

file = '../data/2607_mc_hmc_perc_meth_test.txt'
d = read_bed(file = file, genome = 'hg19', stranded = TRUE, use.score = TRUE)
annotations = c('cpgs', 'detailed_genes')
i = intersect_annotations(
  regions = d,
  annotations = annotations,
  genome = 'hg19',
  ignore.strand = F)
t = annotate_intersections(
  regions = d,
  intersections = i,
  use.score = T)

# Aggregate the score over
aggregation_over_features = dplyr::summarize(
  dplyr::group_by(t, annot_type, annot_id),
  mean = mean(data_score),
  sd = sd(data_score))
aggregation_over_features[[1]] = factor(aggregation_over_features[[1]],
  levels=c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter','hg19_knownGenes_1to5kb','hg19_knownGenes_promoters','hg19_knownGenes_exons5UTRs','hg19_knownGenes_introns5UTRs','hg19_knownGenes_exonsCDSs','hg19_knownGenes_intronsCDSs','hg19_knownGenes_exons3UTRs','hg19_knownGenes_introns3UTRs'))

plot =
  ggplot(aggregation_over_features, aes(mean)) +
  geom_histogram(binwidth=5, aes(y=..density..)) +
  facet_wrap(~annot_type) +
  theme_bw()
ggsave(file = '../figures/2607_mc_hmc_mean_perc_meth_over_annots.png',
  plot = plot, width = 12, height = 12, dpi = 300)

plot =
  ggplot(aggregation_over_features, aes(sd)) +
  geom_histogram(binwidth=5, aes(y=..density..)) +
  facet_wrap(~annot_type) +
  theme_bw()
ggsave(file = '../figures/2607_mc_hmc_sd_perc_meth_over_annots.png',
  plot = plot, width = 12, height = 12, dpi = 300)

################################################################################
# Example of regions with classification and no straind (HUGE, ~25.4M rows)
# Runs read to vis in ~30m (base)
# Runs read to annotate in 2:55m (dplyr)

file = '../data/IDH2mut_v_NBM_class_comp_trim.bed'
d = read_bed(file = file, genome = 'hg19', stranded = FALSE, use.score = FALSE)
annotations = c('cpgs', 'detailed_genes')
i = intersect_annotations(
  regions = d,
  annotations = annotations,
  genome = 'hg19',
  ignore.strand = T)
t = annotate_intersections(
  regions = d,
  intersections = i,
  use.score = F)

plot =
  ggplot(t, aes(data_name, fill=annot_type)) +
  geom_bar()
ggsave(file = '../figures/IDH2mut_v_NBM_class_comp_by_type.png',
  plot = plot, width = 12, height = 12, dpi = 300)

plot =
  ggplot(subset(t, !(data_name %in% c('unclassifiable','noDM'))), aes(data_name, fill=annot_type)) +
  geom_bar()
ggsave(file = '../figures/IDH2mut_v_NBM_class_comp_by_type_no_noDM_uncl.png',
  plot = plot, width = 12, height = 6, dpi = 300)
