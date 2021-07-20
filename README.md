# MapInsights
MapInsights is an efficient toolkit that perform quality control (QC) analysis of next generation sequencing data. Four different modules namely <b><i>`bamqc`</i></b>, <b><i>`genedepth`</i></b>, <b><i>`siteinfo`</i></b> & <b><i>`jumpreads`</i></b> are designed and developed as a part of this toolkit.

<b><i>bamqc</i></b> in an efficient analytical module that performs QC analysis of alignment files. Along with standard QC metrics, it computes extended sets of logs and plots on base mismatches with respect to reference nucleotide in an overall and context specific manner like combination of reads and strands. A coordinate sorted BAM file and reference file are the primary input to the module. Logs, images and a static report in HTML contains summary statistics and plots are the output of <i>bamqc</i>.

<b><i>genedepth</i></b> module calculate exon or region wise depth of coverage. The module takes genomic coordinates in bed format, a sorted BAM file and reference file as input and generate logs and depth plots which are presented as a report in HTML.

<b><i>siteinfo</i></b> module query genomic locus in BAM files and provide comprehensive information about the alignment in query sites such as which nucleotides present in that site, their base quality, strands, read mapping quality, insert-size, read-group etc.
