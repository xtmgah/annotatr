devtools::load_all()
library(ggplot2)

################################################################################
# Section 3.2

hpv_p = '../data/hpv+_mc_hmc_perc_meth.txt'
hpv_n = '../data/hpv-_mc_hmc_perc_meth.txt'

rp = read_bed(hpv_p, genome = 'hg19', stranded = T, use.score = T)
rn = read_bed(hpv_n, genome = 'hg19', stranded = T, use.score = T)

ap = annotate_regions(rp, annotations = c('cpgs','detailed_genes'), genome = 'hg19', ignore.strand = TRUE, use.score = TRUE)
an = annotate_regions(rn, annotations = c('cpgs','detailed_genes'), genome = 'hg19', ignore.strand = TRUE, use.score = TRUE)

sp = summarize_score(ap)
sn = summarize_score(an)

cpgs_order = c(
  'hg19_cpg_islands',
  'hg19_cpg_shores',
  'hg19_cpg_shelves',
  'hg19_cpg_inter')
vp_cpg = visualize_score(sp, cpgs_order)
vn_cpg = visualize_score(sn, cpgs_order)

ggsave(filename='../paper/hpv+_score_over_cpgs.png',width=6, height=6, plot=vp_cpg, dpi=300)
ggsave(filename='../paper/hpv-_score_over_cpgs.png',width=6, height=6, plot=vn_cpg, dpi=300)

genes_order = c(
  'hg19_knownGenes_1to5kb',
  'hg19_knownGenes_promoters',
  'hg19_knownGenes_exons5UTRs',
  'hg19_knownGenes_introns5UTRs',
  'hg19_knownGenes_exonsCDSs',
  'hg19_knownGenes_intronsCDSs',
  'hg19_knownGenes_exons3UTRs',
  'hg19_knownGenes_introns3UTRs')
vp_genes = visualize_score(sp, genes_order)
vn_genes = visualize_score(sn, genes_order)

ggsave(filename='../paper/hpv+_score_over_detailed_genes.png',width=12, height=12, plot=vp_genes, dpi=300)
ggsave(filename='../paper/hpv-_score_over_detailed_genes.png',width=12, height=12, plot=vn_genes, dpi=300)

################################################################################
# Section 3.3

dm = '../data/IDH2mut_v_NBM_DM.txt'

rdm = read_bed(dm, genome = 'hg19', stranded = F, use.score = T)
adm = annotate_regions(rdm, annotations = c('cpgs','detailed_genes'), genome = 'hg19', ignore.strand = TRUE, use.score = TRUE)
sdm = summarize_name(adm)

x_order = c('DMup','DMdown')
fill_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter')
vdm_cpgs = visualize_name(
  summarized_names = sdm,
  x = 'data_name',
  fill = 'annot_type',
  x_order = x_order,
  fill_order = fill_order,
  position = 'stack'
  )

ggsave(filename='../paper/DM_in_cpgs.png',width=6, height=6, plot=vdm_cpgs, dpi=300)

vdm_cpgs_fill = visualize_name(
  summarized_names = sdm,
  x = 'data_name',
  fill = 'annot_type',
  x_order = x_order,
  fill_order = fill_order,
  position = 'fill'
  )

ggsave(filename='../paper/DM_in_cpgs_fill.png',width=6, height=6, plot=vdm_cpgs_fill, dpi=300)


fill_order = c(
  'hg19_knownGenes_1to5kb',
  'hg19_knownGenes_promoters',
  'hg19_knownGenes_exons5UTRs',
  'hg19_knownGenes_introns5UTRs',
  'hg19_knownGenes_exonsCDSs',
  'hg19_knownGenes_intronsCDSs',
  'hg19_knownGenes_exons3UTRs',
  'hg19_knownGenes_introns3UTRs')
x_order = c('DMup','DMdown')
vdm_genes = visualize_name(
  summarized_names = sdm,
  x = 'data_name',
  fill = 'annot_type',
  x_order = x_order,
  fill_order = fill_order,
  position = 'stack'
  )

ggsave(filename='../paper/DM_in_genes.png',width=6, height=6, plot=vdm_genes, dpi=300)

vdm_genes_fill = visualize_name(
  summarized_names = sdm,
  x = 'data_name',
  fill = 'annot_type',
  x_order = x_order,
  fill_order = fill_order,
  position = 'fill'
  )

ggsave(filename='../paper/DM_in_genes_fill.png',width=6, height=6, plot=vdm_genes_fill, dpi=300)

################################################################################
# Section 3.4

cl = '../data/IDH2mut_v_NBM_class_comp_trim.bed'

rcl = read_bed(cl, genome = 'hg19', stranded = F, use.score = F)
acl = annotate_regions(rcl, annotations = c('cpgs','detailed_genes'), genome = 'hg19', ignore.strand = TRUE, use.score = FALSE)
scl = summarize_name(acl)

x_order = c(
  'hyper5mC_5hmC',
  'hypo5mC_5hmC',
  'hyper5hmC',
  'hypo5hmC',
  'hyper5mC',
  'hypo5mC',
  'hyper5mC_hypo5hmC',
  'hypo5mC_hyper5hmC')
fill_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter')
vcl_cpgs = visualize_name(
  summarized_names = scl,
  x = 'data_name',
  fill = 'annot_type',
  x_order = x_order,
  fill_order = fill_order,
  position = 'stack'
  )

ggsave(filename='../paper/classes_in_cpgs.png',width=6, height=6, plot=vcl_cpgs, dpi=300)

vcl_cpgs_fill = visualize_name(
  summarized_names = scl,
  x = 'data_name',
  fill = 'annot_type',
  x_order = x_order,
  fill_order = fill_order,
  position = 'fill'
  )

ggsave(filename='../paper/classes_in_cpgs_fill.png',width=6, height=6, plot=vcl_cpgs_fill, dpi=300)


fill_order = c(
  'hg19_knownGenes_1to5kb',
  'hg19_knownGenes_promoters',
  'hg19_knownGenes_exons5UTRs',
  'hg19_knownGenes_introns5UTRs',
  'hg19_knownGenes_exonsCDSs',
  'hg19_knownGenes_intronsCDSs',
  'hg19_knownGenes_exons3UTRs',
  'hg19_knownGenes_introns3UTRs')
x_order = c(
  'hyper5mC_5hmC',
  'hypo5mC_5hmC',
  'hyper5hmC',
  'hypo5hmC',
  'hyper5mC',
  'hypo5mC',
  'hyper5mC_hypo5hmC',
  'hypo5mC_hyper5hmC')
scl_genes = visualize_name(
  summarized_names = scl,
  x = 'data_name',
  fill = 'annot_type',
  x_order = x_order,
  fill_order = fill_order,
  position = 'stack'
  )

ggsave(filename='../paper/classes_in_genes.png',width=6, height=6, plot=scl_genes, dpi=300)

scl_genes_fill = visualize_name(
  summarized_names = scl,
  x = 'data_name',
  fill = 'annot_type',
  x_order = x_order,
  fill_order = fill_order,
  position = 'fill'
  )

ggsave(filename='../paper/classes_in_genes_fill.png',width=6, height=6, plot=scl_genes_fill, dpi=300)
