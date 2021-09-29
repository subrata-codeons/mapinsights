// Need to do

 #include <stdio.h>
 #include <stdlib.h>
 #include <time.h>
 #include "bam.h"
 


 #define MAXREADLENGTH 250
 #define MAX_ISIZE 600
 #define BAM_FSUPPLEMENTARY 2048

 
typedef struct {
    char *out_file, *bam_file;
    int min_isize, max_ovrlap;
} paramtr;


void jumpreads_help()
 {
   fprintf(stderr, "\n");
   fprintf(stderr, "Usage: mapinsights jumpreads -o <output-file> -i <aligned.bam>\n\n");
   fprintf(stderr, "Options:\n");
   fprintf(stderr, "	-s  minimum insertsize	[1k]\n");
   fprintf(stderr, "	-i  input file	Alignment file (BAM)\n");
   fprintf(stderr, "	-o  output alignment file (BAM)\n");
   fprintf(stderr, "	-c  maximum overlap between outward pair(range : [0 to readlength-1])	[90 bases]\n");
   fprintf(stderr, "\n");
   //return 1;
 }


 int main_jumpreads(int argc,char *argv[])
 {
   bamFile fp, ofp;
   bam_header_t *xheader;
   bam1_t *b;
   int iSize, x=0, n=0, rlen=0, olap=0, f;
   int Isize_ge1K=0, rr=0, ff=0, oo=0, map_difchr=0;
   char *out="Output.bam";
   
   paramtr *prmtr = calloc(1, sizeof(paramtr));
   
   while ((n = getopt(argc, argv, "i:o:s:c:h")) >= 0) {
		switch (n) {
			 case 'i': prmtr->bam_file = optarg; break; // Input bam file
			 case 'o': prmtr->out_file = optarg; break; // output folder path
                        case 's': prmtr->min_isize = atoi(optarg); break;
                        case 'c': prmtr->max_ovrlap = atoi(optarg); break;
                        case 'h': jumpreads_help(); return 1; break;
                        }
	}
	if (argc == 1 || optind-1 == argc) {jumpreads_help(); return 1;}
   
   if(access(prmtr->bam_file, F_OK) == -1) {printf("Alignment(.bam) file does not exist.\n"); return 1;}
   
   if(!prmtr->out_file) {prmtr->out_file=out;}
   
   if(prmtr->min_isize) {if(prmtr->min_isize <= 0) {printf("Please provide valid insert size.\n"); return 1;}}
   else {prmtr->min_isize = 1000;}
   
   if(prmtr->max_ovrlap) {if(prmtr->max_ovrlap < 0) {printf("Please provide valid maximum overlap (bp).\n"); return 1;}}
   else {prmtr->max_ovrlap = 90;}
   
   
   fp = bam_open(prmtr->bam_file, "r");
   ofp = bam_open(prmtr->out_file, "w");
   //printf("I am here1\n");
   if(fp==NULL) { printf("Unable to open input BAM file\n"); return 1;}
   if(ofp==NULL) { printf("Unable to open output BAM file\n"); return 1;}
   xheader=bam_header_init();
   xheader=bam_header_read(fp);
   b=bam_init1();
   bam_header_write(ofp, xheader); 
      while( (x=bam_read1(fp,b)) >= 0 )
      {
        rlen = b->core.l_qseq;
        if (b->core.flag & BAM_FUNMAP || b->core.flag & BAM_FSUPPLEMENTARY || b->core.flag & BAM_FSECONDARY) {continue;}
        f=0;
        
	if((b->core.flag & BAM_FREVERSE) && (b->core.flag & BAM_FMREVERSE)) {f=1; rr++;}
	else if(!(b->core.flag & BAM_FREVERSE) && !(b->core.flag & BAM_FMREVERSE)) {f=1; ff++;}
	olap = rlen - prmtr->max_ovrlap;
	if(olap > 0) {
	if(!(b->core.flag & BAM_FREVERSE) && (b->core.flag & BAM_FMREVERSE) && (b->core.mpos - b->core.pos) < ((-1)*olap) && (b->core.tid == b->core.mtid)) {f=1; oo++;}
	else if((b->core.flag & BAM_FREVERSE) && !(b->core.flag & BAM_FMREVERSE) && (b->core.pos - b->core.mpos) < ((-1)*olap) && (b->core.tid == b->core.mtid)) {f=1; oo++;}
	}
	if (b->core.tid != b->core.mtid) {f=1; map_difchr++;}
	
	if(b->core.isize < 0) {iSize=(b->core.isize * (-1));} else {iSize=b->core.isize;}
	if(iSize >= prmtr->min_isize) {Isize_ge1K++; f=1;}
        if(f) { bam_write1(ofp, b); }
      }
      printf("Pair orientation::\n");
      printf("Mapped-pair forward = %d\n", ff);
      printf("Mapped-pair reverse = %d\n", rr);
      printf("Mapped-pair outward = %d\n", oo);
      printf("Insertsize>=%.2fK = %d\n",((float)prmtr->min_isize/1000.00), Isize_ge1K);
      printf("Mapped-to-diffChr = %d\n", map_difchr);
bam_close(fp);
bam_close(ofp);
bam_index_build(prmtr->out_file);
return 0;
}

