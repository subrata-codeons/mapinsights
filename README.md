# MapInsights



MapInsights is an efficient toolkit that perform quality control (QC) analysis of next generation sequencing data. Four different modules namely <b><i>`bamqc`</i></b>, <b><i>`genedepth`</i></b>, <b><i>`siteinfo`</i></b> & <b><i>`jumpreads`</i></b> are designed and developed as a part of this toolkit. Tow other modules <b><i>`multisample-bamqc`</i></b>, <b><i>`batchplot-bamqc`</i></b> also included in our toolkit that work in relation with <b><i>`bamqc`</i></b> module.

<b><i>bamqc</i></b> is an efficient analytical module that performs QC analysis of alignment files. Along with standard QC metrics, it computes extended sets of logs and plots on base mismatches with respect to reference nucleotide in an overall and context specific manner like combination of reads and strands. A coordinate sorted BAM file and reference file are the primary input to the module. Logs, images and a static report in HTML contains summary statistics and plots are the output of <i>bamqc</i>.

<b><i>multisample-bamqc</i></b> perform combined analysis on multisample <b><i>bamqc</i></b> results. The module facilitates comparison analysis across multiple samples that allow evaluation of consistency among samples and outlier detection based on different features.Given a set of multisample <b><i>bamqc</i></b> results the module perform hierarchical cluster analysis to estimate the variability. PDF file containging plots (features and tree) and other log files are the output of <i>multisample-bamqc</i>.


<b><i>batchplot-bamqc</i></b> generate batch-plots based on different features. Given multisample <b><i>bamqc</i></b> results in batch-wise manner it output batch-plot of two given group and put them in side-by-side fashion for better comparison. 

<b><i>genedepth</i></b> module calculate exon or region wise depth of coverage. The module takes genomic coordinates in bed format, a sorted BAM file and reference file as input and generate logs and depth plots which are presented as a report in HTML.

<b><i>siteinfo</i></b> module query genomic locus in BAM files and provide comprehensive information about the alignment in query sites such as which nucleotides present in that site, their base quality, strands, read mapping quality, insert-size, read-group etc. Referecne fasta, BAM files and genomic coordinate(s) are the main inputs of this module. Genomic coordinate can be passed as argument or a text file contains a list of coordinates can be used for batch query.

<b><i>jumpreads</i></b> module extract reads with atypical alignment properties such as extra long inserts, mate map to different contigs and exception in read orientation etc. jumpreads acts on a coordinate sorted BAM file which is the primary input of the module and generate a bam file as output.

# Installation
### **Requirements:**

<b>- R</b>  <i>` -ggplot2, -ggdendro, -gridExtra `</i>

### **Getting started:**
```
git clone https://github.com/subrata-codeons/mapinsights
cd mapinsights; make
```
OR

download the comprass file and
```
unzip mapinsights-master
cd mapinsights-master; make
```
# Usage
### **mapinsights**
```
Program: mapinsights
Version: 1.0

Usage:   mapinsights <command> [options]

Command: bamqc        QC of alignment file
	 multisample-bamqc	Multi-sample QC analysis based on bamqc results
         batchplot-bamqc	Plotting multi-sample bamqc results
         genedepth    Estimate exon-wise bed coverage
         siteinfo     Details information about genomic site(s)
         jumpreads    Extract reads with jump alignment
```

### **mapinsights-bamqc**
```
Usage:  mapinsights bamqc -r <ref.fa> -o <output-folder-path> -i <aligned.bam>

Options:
	-r   ref.fa                        reference fasta
	-b   bed file                      regions (BED) [null]
	-i   input file                    Alignment file (BAM)
	-o   output folder                 [./]
	-x   exclude read groups
             listed in FILE, one per line  [null]
	-h   help
```
### **mapinsights-multisample-bamqc**
```
Usage:  mapinsights multisample-bamqc -i <folder_path.list>

Options:
	-r   folder_path.list
	-h   help

About input file structure:
folder_path.list file containing a list of bamqc output folder path, one per line.

There are two tab separated columns in the input file (folder_path.list), 
first column should be sample_id and the second column should be the path to 
bamqc output folder. Example structure of the file is given below:

 [sample1	path_to_bamqc_output_folder_of_sample1]
 [sample2	path_to_bamqc_output_folder_of_sample2]
 [  .	       .                                      ]
 [  .	       .                                      ]
 [sampleN	path_to_bamqc_output_folder_of_sampleN]
```
### **mapinsights-batchplot-bamqc**
```
Usage:  mapinsights multisample-bamqc -1 <bamqc_folder_path_for_batch1.list> -2 <bamqc_folder_path_for_batch2.list>

Options:
	-1   bamqc_folder_path_for_batch1.list
	-2   bamqc_folder_path_for_batch2.list
	-h   help
About input file structure:
folder_path_for_batch1/2.list file containing a list of bamqc output folder path, one per line.

There are two tab separated columns in the input file (folder_path.list), 
first column should be sample_id and the second column should be the path to 
bamqc output folder. Example structure of the file is given below:

 [sample1	path_to_bamqc_output_folder_of_sample1]
 [sample2	path_to_bamqc_output_folder_of_sample2]
 [  .	       .                                      ]
 [  .	       .                                      ]
 [sampleN	path_to_bamqc_output_folder_of_sampleN]
```
### **mapinsights-genedepth**
```
Usage:  mapinsights genedepth -r <ref.fa> -o <output-folder> -i <aligned.bam> -b <gene.bed>

Options:
	-r  ref.fa              reference fasta
	-b  gene bed file	regions (BED)
	-i  input file          Alignment file (BAM)
	-o  output folder	[./] 
	-h  help
```
### **mapinsights-siteinfo**
```
Usage: mapinsights siteinfo -r <ref.fa> -i <aligned.bam> 

Options:
	-r  ref.fa                	reference fasta
	-i  input file            	Alignment file (BAM)
	-o  output file 
	-s  query coordinate (single)   [chr#:position-position]
	-p  query coordinate (batch)
	    listed in FILE, one per line     [example :chr#	pos]
	-h  help

```
### **mapinsights-jumpreads**
```
Usage: mapinsights jumpreads -o <output-file> -i <aligned.bam>

Options:
	-s  minimum insertsize		[1k]
	-i  input file			Alignment file (BAM)
	-o  output 			Alignment file (BAM)
	-c  maximum overlap between outward pair(range : [0 to readlength-1])	[90 bases]
	

```
[Mapinsights-bamqc module report (example)](https://subrata-codeons.github.io/Demo/Bamqc.html#GC%20content%20distribution)

[Mapinsights-genedepth module report (example)](https://subrata-codeons.github.io/Demo/Gene_depth.html#GC%20content%20distribution)


