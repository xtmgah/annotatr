## Description of `data-raw/`

This is where data downloaded from data sources goes. We currently support human genome versions hg38, hg19, and mouse genome versions mm10 and mm9.

### `*_cpg_islands.txt.gz`
CpG island tracks are from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables).

### `*_knownGenes.txt.gz`
UCSC knownGene tracks are from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables).

### `*.chrom.sizes.txt`
Chromosome sizes are from [UCSC Genome Browser](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&chromInfoPage=). In the URL, `db=hg38` is replaced with the desired genome, and the `*.chrom.sizes.txt` is downloaded.
