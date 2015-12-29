# This code constructs files relating to introns and exons
# NOTE: This script is run in data-raw/build_data.sh
# DO NOT RUN ALONE

kg_files = list.files(pattern='coding', full.names=T)

for(file in kg_files) {
  message(sprintf('Processing exons/introns for %s', file))

  # Get genome and create temporary file names
    genome = unlist(strsplit(gsub('tmp_', '', basename(file)), '_'))[1]
    exons_file = gsub('coding', 'exons', file)
    introns_file = gsub('coding', 'introns', file)
    firstexons_file = gsub('coding', 'firstexons', file)
    firstintrons_file = gsub('coding', 'firstintrons', file)

  # Get tmp data
    tmp_data = read.table(file, sep='\t', header=F, stringsAsFactors=F)
    colnames(tmp_data) = c(
      'name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd',
      'exonCount','exonStarts','exonEnds','proteinID','alignID')

  # Go through data and build intron and exon lists
    for(i in 1:nrow(tmp_data)) {
      if(i %% 5000 == 0) message(sprintf('On transcript %s of %s', i, nrow(tmp_data)))

      tmp = tmp_data[i,]

      tmp_exonstarts = unlist(strsplit(tmp$exonStarts, ','))
      tmp_exonends = unlist(strsplit(tmp$exonEnds, ','))

      if(tmp$strand == '+') {
        firstexons = data.frame(
          chrom = tmp$chrom,
          start = tmp_exonstarts[1],
          end = tmp_exonends[1],
          name = tmp$name,
          strand = tmp$strand,
          stringsAsFactors=F
        )
      } else if (tmp$strand == '-') {
        firstexons = data.frame(
          chrom = tmp$chrom,
          start = tmp_exonstarts[length(tmp_exonstarts)],
          end = tmp_exonends[length(tmp_exonends)],
          name = tmp$name,
          strand = tmp$strand,
          stringsAsFactors=F
        )
      }
      write.table(firstexons, file=firstexons_file, sep='\t', quote=F, append=T, col.names=F, row.names=F)

      exons = data.frame(
        chrom = tmp$chrom,
        start = tmp_exonstarts,
        end = tmp_exonends,
        name = tmp$name,
        strand = tmp$strand,
        stringsAsFactors=F
      )
      write.table(exons, file=exons_file, sep='\t', quote=F, append=T, col.names=F, row.names=F)

      if(nrow(exons) > 1) {
        if(tmp$strand == '+') {
          firstintrons = data.frame(
            chrom = tmp$chrom,
            start = as.integer(exons$end[1]),
            end = as.integer(exons$start[2]),
            name = tmp$name,
            strand = tmp$strand,
            stringsAsFactors=F
          )
        } else if (tmp$strand == '-') {
          firstintrons = data.frame(
            chrom = tmp$chrom,
            start = as.integer(exons$end[length(exons$end)-1]),
            end = as.integer(exons$start[length(exons$start)]),
            name = tmp$name,
            strand = tmp$strand,
            stringsAsFactors=F
          )
        }
        write.table(firstintrons, file=firstintrons_file, sep='\t', quote=F, append=T, col.names=F, row.names=F)

        introns = data.frame(
          chrom = tmp$chrom,
          start = as.integer(exons$end[-length(exons$end)]),
          end = as.integer(exons$start[-1]),
          name = tmp$name,
          strand = tmp$strand,
          stringsAsFactors=F
        )
        write.table(introns, file=introns_file, sep='\t', quote=F, append=T, col.names=F, row.names=F)
      }
    }
}
