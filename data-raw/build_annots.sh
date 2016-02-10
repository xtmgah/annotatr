#!/bin/bash
set -e
set -u
set -o pipefail

# Usage: bash build_data.sh builds tmp files and runs build_data.R to Construct
# GenomicRanges objects in data/

################################################################
# hg19 Enhancers
################################################################

echo Processing hg19 Fantom5 Permissive Enhancers
awk -v OFS='\t' '{print $1, $2, $3, "enhancer:"NR}' hg19_enhancers_fantom.bed \
| sort -T . -k1,1 -k2,2n \
> renhtmp_hg19_enhancers_fantom.txt

# Process the CpG annotations
for genome in {'hg19','hg38','mm9','mm10'}
do
  ################################################################
  # CpG annotations
  ################################################################
  echo Processing CpG annotations: ${genome}

  # Process CpG islands
  echo Processing ${genome} CpG islands
  awk -v OFS='\t' 'NR > 1 {print $2, $3, $4}' <(gunzip -c ${genome}_cpg_islands.txt.gz) \
  | sort -T . -k1,1 -k2,2n \
  | awk -v OFS='\t' '{print $1, $2, $3, "island:"NR}' \
  > rcpgtmp_${genome}_cpg_islands.txt

  # Process CpG shores
  echo Processing ${genome} CpG shores
  bedtools subtract \
    -a <(bedtools flank -b 2000 -i rcpgtmp_${genome}_cpg_islands.txt -g ${genome}.chrom.sizes.txt | sort -T . -k1,1 -k2,2n | bedtools merge) \
    -b rcpgtmp_${genome}_cpg_islands.txt \
  | awk -v OFS='\t' '{print $1, $2, $3, "shore:"NR}' \
  > rcpgtmp_${genome}_cpg_shores.txt

  # Process CpG shelves
  echo Processing ${genome} CpG shelves
  bedtools subtract \
    -a <(bedtools flank -b 2000 -i rcpgtmp_${genome}_cpg_shores.txt -g ${genome}.chrom.sizes.txt | sort -T . -k1,1 -k2,2n | bedtools merge) \
    -b <(cat rcpgtmp_${genome}_cpg_islands.txt rcpgtmp_${genome}_cpg_shores.txt | sort -T . -k1,1 -k2,2n) \
  | awk -v OFS='\t' '{print $1, $2, $3, "shelf:"NR}' \
  > rcpgtmp_${genome}_cpg_shelves.txt

  # Process inter CpG annotations
  echo Processing ${genome} interCGI
  bedtools complement \
    -i <(cat rcpgtmp_${genome}_cpg_islands.txt rcpgtmp_${genome}_cpg_shores.txt rcpgtmp_${genome}_cpg_shelves.txt | sort -T . -k1,1 -k2,2n | bedtools merge) \
    -g <(sort -k1,1 -k2,2n ${genome}.chrom.sizes.txt) \
  | awk -v OFS='\t' '$2 != $3 {print $1, $2, $3, "inter:"NR}' \
  > rcpgtmp_${genome}_cpg_inter.txt

  ################################################################
  # knownGenes
  # NOTE: Strandedness makes starts/ends on + strand are correct while
  # starts/ends on - strand are reversed!
  ################################################################

  # Order of columns in knownGenes files
  # name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	proteinID	alignID
  #   1     2     3       4       5      6         7        8         9          10        11       12

  echo Processing knownGenes: ${genome}

  # Pull out + and - strands of coding transcripts
  # Do not allow transcripts to start before 5000Kb
  echo Processing ${genome} coding positive strand transcripts
  awk -v OFS='\t' 'NR > 1 && $6 != $7 && $3 == "+" && $4 > 5000 {print $0}' <(gunzip -c ${genome}_knownGenes.txt.gz) \
  | sort -T . -k1,1 -k2,2n \
  > tmp_${genome}_coding_pos_strand_knownGenes.txt

  echo Processing ${genome} coding negative strand transcripts
  awk -v OFS='\t' 'NR > 1 && $6 != $7 && $3 == "-" && $4 > 5000 {print $0}' <(gunzip -c ${genome}_knownGenes.txt.gz) \
  | sort -T . -k1,1 -k2,2n \
  > tmp_${genome}_coding_neg_strand_knownGenes.txt

  # For positive strand features, starts are starts and ends are ends
  # For negative strand features, starts are ends and ends are starts
  for strand in {'pos','neg'}
  do
    echo Processing ${strand} strand
    if [ $strand == 'pos' ]
    then
      # Pull out 5'UTRs, when they exist, from knownGenes
      # chrom, start, end, name, strand
      echo Processing ${genome} ${strand} 5-prime UTRs
      awk -v OFS='\t' '$4 != $6 {print $2, $4, $6, $1, $3}' tmp_${genome}_coding_${strand}_strand_knownGenes.txt \
      | sort -T . -k1,1 -k2,2n \
      > tmp_${genome}_5UTRs_${strand}_strand_knownGenes.txt

      # Pull out 3'UTRs, when they exist, from knownGenes
      # chrom, start, end, name, strand
      echo Processing ${genome} ${strand} 3-prime UTRs
      awk -v OFS='\t' '$5 != $7 {print $2, $7, $5, $1, $3}' tmp_${genome}_coding_${strand}_strand_knownGenes.txt \
      | sort -T . -k1,1 -k2,2n \
      > tmp_${genome}_3UTRs_${strand}_strand_knownGenes.txt

      # Define promoter as <= 1Kb upstream of TSS
      # chrom, start, end, name, strand
      echo Processing ${genome} ${strand} promoters
      bedtools flank \
        -l 1000 -r 0 \
        -i <(awk -v OFS='\t' '{print $2, $4, $4, $1, $3}' tmp_${genome}_coding_${strand}_strand_knownGenes.txt) \
        -g ${genome}.chrom.sizes.txt \
      | sort -T . -k1,1 -k2,2n \
      | awk -v OFS='\t' '$2 < $3 {print $0}' \
      > tmp_${genome}_promoters_${strand}_strand_knownGenes.txt

      # Pull out region > 1Kb and <= 5Kb upstream of TSS
      # chrom, start, end, name, strand
      echo Processing ${genome} ${strand} 1to5kb upstream of TSS
      bedtools flank \
        -l 4000 -r 0 \
        -i <(awk -v OFS='\t' '{print $1, $2, $2, $4, $5}' tmp_${genome}_promoters_${strand}_strand_knownGenes.txt) \
        -g ${genome}.chrom.sizes.txt \
      | sort -T . -k1,1 -k2,2n \
      | awk -v OFS='\t' '$2 < $3 {print $0}' \
      > tmp_${genome}_1to5kb_${strand}_strand_knownGenes.txt
    else
      # Pull out 5'UTRs, when they exist, from knownGenes
      # chrom, start, end, name, strand
      echo Processing ${genome} ${strand} 5-prime UTRs
      awk -v OFS='\t' '$5 != $7 {print $2, $7, $5, $1, $3}' tmp_${genome}_coding_${strand}_strand_knownGenes.txt \
      | sort -T . -k1,1 -k2,2n \
      > tmp_${genome}_5UTRs_${strand}_strand_knownGenes.txt

      # Pull out 3'UTRs, when they exist, from knownGenes
      # chrom, start, end, name, strand
      echo Processing ${genome} ${strand} 3-prime UTRs
      awk -v OFS='\t' '$4 != $6 {print $2, $4, $6, $1, $3}' tmp_${genome}_coding_${strand}_strand_knownGenes.txt \
      | sort -T . -k1,1 -k2,2n \
      > tmp_${genome}_3UTRs_${strand}_strand_knownGenes.txt

      # Define promoter as <= 1Kb upstream of TSS
      # chrom, start, end, name, strand
      echo Processing ${genome} ${strand} promoters
      bedtools flank \
        -r 1000 -l 0 \
        -i <(awk -v OFS='\t' '{print $2, $5, $5, $1, $3}' tmp_${genome}_coding_${strand}_strand_knownGenes.txt) \
        -g ${genome}.chrom.sizes.txt \
      | sort -T . -k1,1 -k2,2n \
      | awk -v OFS='\t' '$2 < $3 {print $0}' \
      > tmp_${genome}_promoters_${strand}_strand_knownGenes.txt

      # Pull out region > 1Kb and <= 5Kb upstream of TSS
      # chrom, start, end, name, strand
      echo Processing ${genome} ${strand} 1to5kb upstream of TSS
      bedtools flank \
        -r 4000 -l 0 \
        -i <(awk -v OFS='\t' '{print $1, $3, $3, $4, $5}' tmp_${genome}_promoters_${strand}_strand_knownGenes.txt) \
        -g ${genome}.chrom.sizes.txt \
      | sort -T . -k1,1 -k2,2n \
      | awk -v OFS='\t' '$2 < $3 {print $0}' \
      > tmp_${genome}_1to5kb_${strand}_strand_knownGenes.txt
    fi

    # Pull out CDSs
    echo Processing ${genome} ${strand} CDSs
    awk -v OFS='\t' '{print $2, $6, $7, $1, $3}' tmp_${genome}_coding_${strand}_strand_knownGenes.txt \
    > tmp_${genome}_CDSs_${strand}_strand_knownGenes.txt
  done

done

# Process exons and introns
Rscript build_annots_exint.R

for genome in {'hg19','hg38','mm9','mm10'}
do
  for strand in {'pos','neg'}
  do
    echo Sorting ${genome} introns/exons on ${strand} strand
    sort -T . -k1,1 -k2,2n tmp_${genome}_exons_${strand}_strand_knownGenes.txt > tmp_${genome}_exons_${strand}_strand_sorted.txt
    mv tmp_${genome}_exons_${strand}_strand_sorted.txt tmp_${genome}_exons_${strand}_strand_knownGenes.txt

    sort -T . -k1,1 -k2,2n tmp_${genome}_introns_${strand}_strand_knownGenes.txt > tmp_${genome}_introns_${strand}_strand_sorted.txt
    mv tmp_${genome}_introns_${strand}_strand_sorted.txt tmp_${genome}_introns_${strand}_strand_knownGenes.txt

    echo Sorting ${genome} first introns / first exons on ${strand} strand
    sort -T . -k1,1 -k2,2n tmp_${genome}_firstexons_${strand}_strand_knownGenes.txt > tmp_${genome}_firstexons_${strand}_strand_sorted.txt
    mv tmp_${genome}_firstexons_${strand}_strand_sorted.txt tmp_${genome}_firstexons_${strand}_strand_knownGenes.txt

    sort -T . -k1,1 -k2,2n tmp_${genome}_firstintrons_${strand}_strand_knownGenes.txt > tmp_${genome}_firstintrons_${strand}_strand_sorted.txt
    mv tmp_${genome}_firstintrons_${strand}_strand_sorted.txt tmp_${genome}_firstintrons_${strand}_strand_knownGenes.txt
  done
done

# More finely process exons and introns
for genome in {'hg19','hg38','mm9','mm10'}
do
  for strand in {'pos','neg'}
  do
    echo Processing ${genome} ${strand}
    for annot in {'5UTRs','3UTRs','CDSs'}
    do
      ########################################################
      # Capture exons UTRs
      # See figure in meta for clarification of cases
      # exonStart exonEnd utrStart utrEnd
      #    2         3       7        8
      echo Processing ${genome} ${strand} ${annot}-exons
      bedtools intersect -wa -wb \
        -a tmp_${genome}_exons_${strand}_strand_knownGenes.txt \
        -b tmp_${genome}_${annot}_${strand}_strand_knownGenes.txt \
      | awk -v OFS='\t' \
        '$4 == $9 && $2 <= $7 && $3 > $7 && $3 <  $8 {print $1, $7, $3, $4, $5} \
         $4 == $9 && $2 <= $7 && $3 >= $8 {print $1, $7, $8, $4, $5} \
         $4 == $9 && $2 >  $7 && $3 <  $8 {print $1, $2, $3, $4, $5} \
         $4 == $9 && $2 >  $7 && $2 < $8 && $3 >= $8 {print $1, $2, $8, $4, $5}'\
      | sort -T . -k1,1 -k2,2n \
      > tmp_${genome}_exons${annot}_${strand}_strand_knownGenes.txt

      # Capture introns in UTRs
      # By definition, the cdsStart/cdsEnd cannot be in an intron
      # nor can the txStart/txEnd be in an intron.
      # intronStart intronEnd utrStart utrEnd
      #      2          3        7       8
      echo Processing ${genome} ${strand} ${annot}-introns
      bedtools intersect -wa -wb \
        -a tmp_${genome}_introns_${strand}_strand_knownGenes.txt \
        -b tmp_${genome}_${annot}_${strand}_strand_knownGenes.txt \
      | awk -v OFS='\t' '$4 == $9 {print $1, $2, $3, $4, $5}' \
      | sort -T . -k1,1 -k2,2n \
      > tmp_${genome}_introns${annot}_${strand}_strand_knownGenes.txt
    done
  done
done

# Combine and sort strand annotations
for genome in {'hg19','hg38','mm9','mm10'}
do
  for annot in {'1to5kb','promoters','5UTRs','exons5UTRs','introns5UTRs','exonsCDSs','intronsCDSs','CDSs','firstexons','exons','firstintrons','introns','3UTRs','exons3UTRs','introns3UTRs'}
  do
    echo Combining ${genome} ${annot} across strands
    # Collapse, over all transcripts, features with identical coordinates
    # This is convoluted, but need to preserve strand
    cat \
      <(bedtools groupby \
        -i tmp_${genome}_${annot}_pos_strand_knownGenes.txt \
        -g 1-3 \
        -c 4 \
        -o collapse \
      | awk -v OFS='\t' '{print $1, $2, $3, "+", $4}') \
      <(bedtools groupby \
        -i tmp_${genome}_${annot}_neg_strand_knownGenes.txt \
        -g 1-3 \
        -c 4 \
        -o collapse \
      | awk -v OFS='\t' '{print $1, $2, $3, "-", $4}') \
    | sort -T . -k1,1 -k2,2n \
    | awk -v OFS='\t' '$2 < $3 {print $0}' \
    > rkgtmp_${genome}_knownGenes_${annot}.txt
    rm tmp_${genome}_${annot}_{pos,neg}_strand_knownGenes.txt
  done
done

# Create intergenic annotations as the complement of existing knownGenes
# annotations + gaps
for genome in {'hg19','hg38','mm9','mm10'}
do
  echo Formatting ${genome} gaps
  awk -v OFS='\t' 'NR > 1 { print $2, $3, $4, "*", "gap:"NR }' <(gunzip -c ${genome}_gaps.txt.gz) > tmp_${genome}_gaps.txt

  echo Concatenating and sorting all knownGenes annotations
  cat tmp_${genome}_gaps.txt rkgtmp_${genome}_knownGenes_*.txt \
  | sort -T . -k1,1 -k2,2n \
  > tmp_${genome}_knownGenes_all.txt
  rm tmp_${genome}_gaps.txt

  echo Merging all knownGenes annotations
  bedtools merge -d 5 -i tmp_${genome}_knownGenes_all.txt \
  | sort -k1,1 -k2,2n > tmp_${genome}_knownGenes_merged.txt
  rm tmp_${genome}_knownGenes_all.txt

  echo Building ${genome} intergenic annotations
  bedtools subtract -a <(awk -v OFS='\t' '{ print $1, "0", $2 }' ${genome}.chrom.sizes.txt) -b tmp_${genome}_knownGenes_merged.txt \
  | sort -T . -k1,1 -k2,2n \
  | awk -v OFS='\t' '$2 < $3 { print $1, $2, $3, "*", "intergenic:"NR }' \
  > rkgtmp_${genome}_knownGenes_intergenic.txt
  rm tmp_${genome}_knownGenes_merged.txt
done

# # Remove the coding tmp files
# find . -name '*coding*' | xargs rm
#
# # After pre-processing, run R to convert objects to GRanges
# Rscript build_annots.R
