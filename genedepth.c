
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include "bam.h"
#include "faidx.h"
typedef struct {     
	bamFile fp;      
	bam_iter_t iter; 
	int min_mapQ, min_len; 
} aux_t;



static int read_bam(void *data, bam1_t *b) 
{
	aux_t *aux = (aux_t*)data; 
	int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
	if (!(b->core.flag&BAM_FUNMAP)) {
		if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
		else if (aux->min_len && bam_cigar2qlen(&b->core, bam1_cigar(b)) < aux->min_len) b->core.flag |= BAM_FUNMAP;
	}
	return ret;
}

 
 void genecov_draw_plots(char *outdir, char *fl, int exn_cnt)
 {
   FILE *gfp1;
   char *outfn=0;
   outfn = calloc(strlen(outdir) + 500, 1);
   sprintf(outfn, "%s/plotsmake.r", outdir);
   
   gfp1=fopen(outfn,"w");
   fprintf(gfp1,"library (\"ggplot2\")\n");

   fprintf(gfp1,"dplot <- read.table(\"%s.txt\",header = TRUE, sep = \"\\t\")\n", fl);
   fprintf(gfp1,"Exone%d <- dplot[,c(\"pos\")]\n",exn_cnt);
   fprintf(gfp1,"png(file=\"%s.png\")\n", fl);
   fprintf(gfp1,"ggplot(data = dplot, aes(x = pos, y = depth, fill=factor(type, levels=c(\"nonref\",\"ref\")))) +  geom_area(stat=\"identity\") + scale_fill_manual(\"type\",values = c(\"brown2\",\"gray30\")) + expand_limits(y = 0) + scale_x_continuous(\"Exon%d\", labels = as.character(Exone%d), breaks = Exone%d) + theme(axis.ticks.x = element_line(size = 2, color=\"orange\") , axis.ticks.length = unit(.3, \"cm\")) + theme(axis.text.x=element_blank())\n",exn_cnt ,exn_cnt ,exn_cnt);
   
   fprintf(gfp1,"dev.off()\n");
   fclose(gfp1);
   *outfn=0;
   sprintf(outfn, "R CMD BATCH %s/plotsmake.r", outdir);
   system(outfn);
   *outfn=0;
   sprintf(outfn, "rm -f %s/plotsmake.r", outdir);
   system(outfn);
   system("rm -f plotsmake.*");
   free(outfn);
 }
 void creathtml(char *outdir, int exn_count, int totesize, float dcov)
 {
   int i, j, k;
   FILE *fp, *ifp;
   char p1[1024], p2[1024];
   char *outfn=0, *infn=0;
   infn = calloc(strlen(outdir) + 500, 1);
   outfn = calloc(strlen(outdir) + 500, 1);
   sprintf(infn, "%s/Gene_depth.log", outdir);
   sprintf(outfn, "%s/Gene_depth.html", outdir);
   fp = fopen(outfn, "w");
   
	fprintf(fp,"<!DOCTYPE html>\n");
	fprintf(fp,"<html lang=\"en\">\n");
	fprintf(fp,"<head>\n");
	fprintf(fp,"<meta charset=\"utf-8\">\n");
	fprintf(fp,"<title>Report(genedepth)</title>\n");
	fprintf(fp,"<style>\n");
	fprintf(fp,"    body{        \n");
	fprintf(fp,"        padding-top: 70px;\n");
	fprintf(fp,"        padding-bottom: 40px;\n");
	fprintf(fp,"        padding-left: 20em;\n");
	fprintf(fp,"        \n");
	fprintf(fp,"    }\n");
	fprintf(fp,"    .container{\n");
	fprintf(fp,"        width: 90%%;\n");
	fprintf(fp,"        margin: 0 auto; \n");
	fprintf(fp,"    }\n");
	fprintf(fp,"    .fixed-header{\n");
	fprintf(fp,"        width: 100%%;\n");
	fprintf(fp,"        position: fixed;        \n");
	fprintf(fp,"        background: #333;\n");
	fprintf(fp,"        padding: 20px 0;\n");
	fprintf(fp,"        color: #fff;\n");
	fprintf(fp,"        left: 0;\n");
	fprintf(fp,"    }\n");
	fprintf(fp,"    .fixed-footer{\n");
	fprintf(fp,"        width: 100%%;\n");
	fprintf(fp,"        position: fixed;        \n");
	fprintf(fp,"        background: #333;\n");
	fprintf(fp,"        padding: 10px 0;\n");
	fprintf(fp,"        color: #fff;\n");
	fprintf(fp,"        left: 0;\n");
	fprintf(fp,"    }\n");
	fprintf(fp,"    .fixed-header{\n");
	fprintf(fp,"        top: 0;\n");
	fprintf(fp,"    }\n");
	fprintf(fp,"    .fixed-footer{\n");
	fprintf(fp,"        bottom: 0;\n");
	fprintf(fp,"    }    \n");
	fprintf(fp,"    /* Some more styles to beautify this example */\n");
	fprintf(fp,"   ul.navbar {\n");
	fprintf(fp,"    position: fixed;\n");
	fprintf(fp,"    top: 4em;\n");
	fprintf(fp,"    left: 1em;\n");
	fprintf(fp,"    width: 19em }\n");
	fprintf(fp,"    \n");
	fprintf(fp,"    .container p{\n");
	fprintf(fp,"        line-height: 200px; /* Create scrollbar to test positioning */\n");
	fprintf(fp,"    }\n");
	fprintf(fp,"</style>\n");
	fprintf(fp,"</head>\n");
	fprintf(fp,"<body>\n");
	fprintf(fp,"    <div class=\"fixed-header\">\n");
        fprintf(fp,"<div class=\"headertitle\">\n");
        fprintf(fp,"<span style=\"color:#FF5733; font-size:125%%; float: left; padding-left: 20px;\"><b>Map</b></span>\n");
        fprintf(fp,"<span style=\"color:#16A085; font-size:125%%; float: left;\"><b>Insights</b></span>\n");
        fprintf(fp,"<mid style=\"color:white; font-size:100%%; float: right; padding-right: 20px;\">Mapinsights genedepth report</mid>\n");
        fprintf(fp,"</div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<ul class=\"navbar\">\n");
	fprintf(fp,"<h2 style=\"color:#797D7F;\">CONTENTS</h2>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Parameters & inputs\" style=\"color:#1C2833;\" >Commands and parameters</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Summary statistics\" style=\"color:#1C2833;\" >Summary statistics</a></li>\n");
	fprintf(fp,"</ul>\n");
	fprintf(fp,"\n");
	fprintf(fp,"        \n");
	fprintf(fp,"    \n");
	fprintf(fp,"  <head>\n");
	fprintf(fp,"<style>\n");
	fprintf(fp,"* {\n");
	fprintf(fp,"  box-sizing: border-box;\n");
	fprintf(fp,"\n");
	fprintf(fp,"}\n");
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,".column {\n");
	fprintf(fp," \n");
	fprintf(fp,"  float: left;\n");
	fprintf(fp,"  width: 50%%;\n");
	fprintf(fp,"  padding: 5px;\n");
	fprintf(fp,"  \n");
	fprintf(fp,"}\n");
	fprintf(fp,"\n");
	fprintf(fp,"/* Clearfix (clear floats) */\n");
	fprintf(fp,".row::after {\n");
	fprintf(fp,"  content: \"\";\n");
	fprintf(fp,"  clear: both;\n");
	fprintf(fp,"  display: table;\n");
	fprintf(fp,"}\n");
	fprintf(fp,"</style>\n");
	fprintf(fp,"</head>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<style>\n");
	fprintf(fp,"\n");
	fprintf(fp,"div.mapsummary{\n");
	fprintf(fp,"\n");
	fprintf(fp,"}\n");
	fprintf(fp,"\n");
	fprintf(fp,"table.mapsummary{\n");
	fprintf(fp,"\n");
	fprintf(fp,"	border: 0px;\n");
	fprintf(fp,"  width: 80%%;\n");
	fprintf(fp,"  vertical-align: top;\n");
	fprintf(fp,"\n");
	fprintf(fp,"}\n");
	fprintf(fp,"\n");
	fprintf(fp,"div.table-mapsummary{\n");
	fprintf(fp,"	margin-left:auto; \n");
	fprintf(fp,"  margin-right:auto; \n");
	fprintf(fp,"	padding-bottom: 20px;\n");
	fprintf(fp,"}\n");
	fprintf(fp,"\n");
	fprintf(fp,"td.column1{\n");
	fprintf(fp,"	width: 50%%;\n");
	fprintf(fp,"}\n");
	fprintf(fp,"\n");
	fprintf(fp,"table.infotable {\n");
	fprintf(fp,"\n");
	fprintf(fp,"	font-size:14px;\n");
	fprintf(fp,"	border-color: #AEB6BF;\n");
	fprintf(fp,"	border-width: 1px;\n");
	fprintf(fp,"	border-collapse: collapse;\n");
	fprintf(fp,"}\n");
	fprintf(fp,"table.infotable th {\n");
	fprintf(fp,"	padding: 8px;\n");
	fprintf(fp,"	border-width: 1px;\n");
	fprintf(fp,"	border-style: solid;\n");
	fprintf(fp,"	border-color: #AEB6BF;\n");
	fprintf(fp,"	background-color:#FBFCFC;\n");
	fprintf(fp,"}\n");
	fprintf(fp,"table.infotable tr {\n");
	fprintf(fp,"	background-color:#FBFCFC;\n");
	fprintf(fp,"}\n");
	fprintf(fp,"table.infotable td {\n");
	fprintf(fp,"	border-width: 1px;\n");
	fprintf(fp,"	padding: 8px;\n");
	fprintf(fp,"	border-style: solid;\n");
	fprintf(fp,"	border-color: #AEB6BF;\n");
	fprintf(fp,"}\n");
	fprintf(fp,"\n");
	fprintf(fp,"</style>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<h2 id=\"Parameters & inputs\" style=\"background-color:#E5E7E9;\"><u>Parameters & inputs</u></h2>\n");
	fprintf(fp,"\n");
   
   
 
      ifp=fopen(infn,"r");
      
      i=0;
      while((p1[i++]=fgetc(ifp)) != '=');
      p1[--i]='\0';
      //ch=fgetc(ifp);
      fgetc(ifp);
      i=0;
      while((p2[i++]=fgetc(ifp)) != '\n');
      p2[--i]='\0';
      
        fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Command line</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>%s</td>\n", p2);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
      
  
      fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Parameters</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
      
      for(j = 0 ; j < 3 ; j++) {
      i=0;
      while((p1[i++]=fgetc(ifp)) != '=');
      p1[--i]='\0';
      //ch=fgetc(ifp);
      fgetc(ifp);
      i=0;
      while((p2[i++]=fgetc(ifp)) != '\n');
      p2[--i]='\0';
      //printf("%s = %s\n",p1,p2);
      
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>%s</td>\n", p1);
	fprintf(fp,"<td class=column2>%s</td>\n", p2);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
      }
      //ch=fgetc(ifp);
      fgetc(ifp);
      fprintf(fp,"</table>\n");
      fprintf(fp,"</div>\n");
      
      fprintf(fp,"<div class=mapsummary>\n");
      fprintf(fp,"<h3>Other information</h3>\n");
      fprintf(fp,"<table class=\"mapsummary infotable\">\n");
      
      for(j = 0 ; j < 3 ; j++) {
      i=0;
      while((p1[i++]=fgetc(ifp)) != '=');
      p1[--i]='\0';
      //ch=fgetc(ifp);
      fgetc(ifp);
      i=0;
      while((p2[i++]=fgetc(ifp)) != '\n');
      p2[--i]='\0';
      
        
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>%s</td>\n", p1);
	fprintf(fp,"<td class=column2>%s</td>\n", p2);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
      }
      fprintf(fp,"</table>\n");
      fprintf(fp,"</div>\n");
      
      fprintf(fp,"<h2 id=\"Summary statistics\" style=\"background-color:#E5E7E9;\"><u>Summary statistics & plots</u></h2>\n");
      
      fprintf(fp,"<div class=mapsummary>\n");
fprintf(fp,"<h3>About bed</h3>\n");
fprintf(fp,"<table class=\"mapsummary infotable\">\n");

fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
fprintf(fp,"<td class=column1>No of exones</td>\n");
fprintf(fp,"<td class=column2>%d</td>\n", exn_count);
fprintf(fp,"</tr>\n");

fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
fprintf(fp,"<td class=column1>Target bed size (bp)</td>\n");
fprintf(fp,"<td class=column2>%d</td>\n", totesize);
fprintf(fp,"</tr>\n");

fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
fprintf(fp,"<td class=column1>Mean depth of coverage</td>\n");
fprintf(fp,"<td class=column2>%.2fX</td>\n", dcov);
fprintf(fp,"</tr>\n");


fprintf(fp,"</table>\n");
      fprintf(fp,"</div>\n");
      
    
      
      
      for(j = 0 ; j < exn_count ; j++) {
      fprintf(fp,"<h3 style=\"background-color:#FEF9E7;\" onmouseover=\"this.style.backgroundColor='#F9E79F';\" onmouseout=\"this.style.backgroundColor='#FEF9E7';\">Exon%d</h3>\n",j+1);
      
      fprintf(fp,"<div class=\"row\">\n");

      fprintf(fp,"<div class=\"column\">\n");
        fprintf(fp,"<div><img width=\"400\" height=\"400\" src=\"plots_gene-depth/Exon%d_depth.png\"></div>\n",j+1);
 fprintf(fp,"</div>\n"); 

 fprintf(fp,"<div class=\"column\">\n");
      
      
      fprintf(fp,"<div class=mapsummary>\n");
	
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
      //ch=fgetc(ifp);
      fgetc(ifp);
      for(k = 0; k < 10 ; k++) {
      
      i=0;
      while((p1[i++]=fgetc(ifp)) != '=');
      p1[--i]='\0';
      //ch=fgetc(ifp);
      fgetc(ifp);
      i=0;
      while((p2[i++]=fgetc(ifp)) != '\n');
      p2[--i]='\0';
      
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>%s</td>\n", p1);
	fprintf(fp,"<td class=column2>%s</td>\n", p2);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
      }
      fprintf(fp,"</table>\n");
      fprintf(fp,"</div>\n");
      fprintf(fp,"</div>\n");
      fprintf(fp,"</div>\n");
      }
   
   fprintf(fp,"\n");
	fprintf(fp,"    <div class=\"fixed-footer\">\n");
	fprintf(fp,"        <div style=\"color:white; font-size:100%%; float: right; padding-right: 20px; class=\"container\">Created by Mapinsights (version 1.0) </div>\n");
	fprintf(fp,"    </div>\n");
	fprintf(fp,"</body>\n");
	fprintf(fp,"</html>\n");
   fclose(ifp);
   fclose(fp);
   free(infn);
   free(outfn);
 }
 
 void Arrange_logs(char *outdir)
 {
   char *outfldr=0;
   outfldr = calloc(strlen(outdir) + 500, 1);
   sprintf(outfldr, "mkdir -p %s/plots_gene-depth", outdir);
   system(outfldr);
   *outfldr=0;  
   sprintf(outfldr, "mv %s/*.png %s/plots_gene-depth", outdir, outdir);
   system(outfldr);

   *outfldr=0;
   sprintf(outfldr, "mkdir -p %s/logs_gene-depth", outdir);
   system(outfldr);
   *outfldr=0;  
   sprintf(outfldr, "mv %s/*.txt %s/logs_gene-depth", outdir, outdir);
   system(outfldr);
   free(outfldr);
 }
 
typedef struct {
 char *ref_file, *bed_file, *bam_file;
 char *out_foldr;
 
} pramtr;

void genedepth_help()
 {
   fprintf(stderr, "\n");
   fprintf(stderr, "Usage: mapinsights genedepth -r <ref.fa> -o <output-folder> -i <aligned.bam> -b <gene.bed>\n\n");
   fprintf(stderr, "Options:\n");
   fprintf(stderr, "	-r  ref.fa      reference fasta\n");
   fprintf(stderr, "	-b  gene bed file	regions (BED)\n");
   fprintf(stderr, "	-i  input file   Alignment file (BAM)\n");
   fprintf(stderr, "	-o  path to output folder	[./] \n");
   fprintf(stderr, "	-h  help\n");
   fprintf(stderr, "\n");
   //return 1;
 }

int main_genedepth(int argc, char *argv[])
{
	int i, n,j, tid, beg, end, pos, prepos, *n_plp;
	const bam_pileup1_t *q,*p;
	bam_header_t *h = 0; 
	aux_t *data;
        faidx_t *fai;
	bam_plp_t mplp;
        int tid0 = -1, ref_tid = -1, ref_len, exn_cnt=1;
        char *ref, ch, tchr[20], tpos_srt[50], tpos_end[50], chrpos[200];
        char *outfn, *opth="./";
        int ref_cnt, var_cnt, dedepth;
        long int ref_size = 0;
        int pct1, pct5, pct10, pct20, pct30;
        int totdepth, esize, gc, totesize = 0, totaldepth = 0;
        FILE *ifp, *ofp, *ofp1;
	pramtr *pmtr = calloc(1, sizeof(pramtr));
	
	n_plp = calloc(1, sizeof(int));
        mplp = calloc(1, sizeof(bam_plp_t));
        while ((n = getopt(argc, argv, "r:b:i:o:h")) >= 0) {
             switch (n) {
                    case 'r': pmtr->ref_file = optarg; break; // reference file
                    case 'b': pmtr->bed_file = optarg; break; // bed file
                    case 'i': pmtr->bam_file = optarg; break; // Input bam file
                    case 'o': pmtr->out_foldr = optarg; break; // output folder path
                    case 'h': genedepth_help(); return 1; break;
                    }
        }
        if (argc == 1 || optind-1 == argc) {genedepth_help(); return 1;}

        if(access(pmtr->ref_file, F_OK) == -1) {printf("Reference file does not exist.\n"); return 1;}

        if(access(pmtr->bam_file, F_OK) == -1) {printf("Alignment(.bam) file does not exist.\n"); return 1;}

        if(access(pmtr->bed_file, F_OK) == -1) {printf("Gene bed file does not exist.\n"); return 1;}

        if(pmtr->out_foldr) {if(access(pmtr->out_foldr, F_OK) == -1) {printf("Output filder does not exist.\n"); return 1;}}
        
        if(!(pmtr->out_foldr)) {pmtr->out_foldr = opth;}
        
        outfn = calloc(strlen(pmtr->out_foldr) + 500, 1);
        
        sprintf(outfn, "%s/Gene_depth.log", pmtr->out_foldr);
        
        ofp1=fopen(outfn,"w");
        
        fprintf(ofp1,"Command = mapinsights");
	for(i=0;i<argc; i++) {fprintf(ofp1," %s", argv[i]);}
	fprintf(ofp1,"\n\nAlignment file = %s\n",pmtr->bam_file);
	fprintf(ofp1,"Reference file = %s\n",pmtr->ref_file);
	fprintf(ofp1,"Bed file = %s\n\n",pmtr->bed_file);
        
	outfn = calloc(strlen(pmtr->out_foldr) + 500, 1);
        
	n=1;
	data = calloc(n, sizeof(void*)); 
	beg = 0; end = 1<<30; tid = -1;  
	        
		bam_header_t *htmp;
               
		data = calloc(1, sizeof(aux_t));
		
		data->fp = bam_open(pmtr->bam_file, "r");
		
		time_t t;
		time(&t);
		fprintf(ofp1,"Analysis date = %s", ctime(&t));
		htmp = bam_header_read(data->fp);        
		for(i=0; i< htmp->n_targets; i++)
               {
                   ref_size += htmp->target_len[i];
        
               }
               fprintf(ofp1,"No of contigs in alignment file = %d\n", htmp->n_targets);
               fprintf(ofp1,"Reference size (bp) = %ld\n\n", ref_size);
                fai = fai_load(pmtr->ref_file);
                tid0 = tid;
        if (tid0 >= 0 && fai) { 
                ref = faidx_fetch_seq(fai, htmp->target_name[tid0], 0, 0x7fffffff, &ref_len);
                ref_tid = tid0;

        } else ref_tid = -1, ref = 0;
                bam_index_t *idx = bam_index_load(pmtr->bam_file);
		h = htmp; 

                ifp=fopen(pmtr->bed_file,"r");
                while((ch=fgetc(ifp))!=EOF)
                {
                   j=0;
                   //printf("I am inside\n");
                   tchr[j]=ch;j++;
                   while(((ch=fgetc(ifp))!='	'))
                   {tchr[j]=ch;j++;}
                   tchr[j]='\0';j=0;
                   //printf("I am inside1\n");
                   while(((ch=fgetc(ifp))!='	'))
                    {tpos_srt[j]=ch;j++;}
                     tpos_srt[j]='\0';j=0;
                   while(((ch=fgetc(ifp))!='	'))
                    {tpos_end[j]=ch;j++;}
                     tpos_end[j]='\0';
                   if(ch != '\n') {while((ch=fgetc(ifp))!='\n');}
                   //printf("%s\n",tchr);
                   //printf("%s\n",tpos);
                   sprintf(chrpos,"%s:%s-%s",tchr,tpos_srt,tpos_end);
                 printf("%s\n",chrpos);
                   
		if (strlen(chrpos) > 0){
                 bam_parse_region(h, chrpos, &tid, &beg, &end);}


		 

		if (tid >= 0) 
		{ 
			
			data->iter = bam_iter_query(idx, tid, beg, end); 
			 
                        
		
	
        if (tid != ref_tid) {
                        free(ref); ref = 0;
                       if (fai) ref = faidx_fetch_seq(fai,  htmp->target_name[tid], 0, 0x7fffffff, &ref_len);
                        ref_tid = tid;}

                }

        pct1=0; pct5=0; pct10=0; pct20=0; pct30=0; totdepth=0; esize=0;
        esize = (end-beg);
        totesize += esize;
        
        gc=0;
        for(i = beg ; i < end; i++) {if(toupper(ref[i]) == 'G' || toupper(ref[i]) == 'C') {gc++;}}
        *outfn=0;
        sprintf(outfn, "%s/Exon%d_depth.txt", pmtr->out_foldr, exn_cnt);
        ofp=fopen(outfn,"w");
        fprintf(ofp,"pos	type	depth\n");
	mplp = bam_plp_init(read_bam, (void*)data); 
	 
	n=1;
        j=0;
	while((q=bam_plp_auto(mplp, &tid, &pos, n_plp)))
        {
          //printf("position = %d\n", pos);
          if (pos < beg || pos >= end) {continue;}
          if(beg < pos && n) {for(i=beg+1; i < pos+1; i++) {fprintf(ofp,"%d	ref	0\n", i); fprintf(ofp,"%d	nonref	0\n", i);}}
          if(!n) {if(pos - prepos > 1) {for(i=prepos+1; i < pos+1; i++) {fprintf(ofp,"%d	ref	0\n", i); fprintf(ofp,"%d	nonref	0\n", i);}}}
          n=0;
          prepos=pos;
	  //printf("%s	%d	depth = %d	", h->target_name[tid], pos+1, *n_plp);
	  ref_cnt=0; var_cnt=0; dedepth=0;
          for (j = 0; j < *n_plp; ++j){ 
          p = &q[j];
          if (p->is_del || p->is_refskip) {dedepth++;}
          else if(bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b),p->qpos)] != toupper(ref[pos])) {var_cnt++;}
          else {ref_cnt++;}
          
          }
          
          if((ref_cnt-dedepth) > 0) 
          {
              fprintf(ofp,"%d	ref	%d\n", pos+1, ref_cnt-dedepth);
              totdepth+=(ref_cnt-dedepth);
              if((ref_cnt-dedepth) >= 1) {pct1++;}
              if((ref_cnt-dedepth) >= 5) {pct5++;}
              if((ref_cnt-dedepth) >= 10) {pct10++;}
              if((ref_cnt-dedepth) >= 20) {pct20++;}
              if((ref_cnt-dedepth) >= 30) {pct30++;}
          } 
          else 
          {
              fprintf(ofp,"%d	ref	0\n", pos+1);
          }
          fprintf(ofp,"%d	nonref	%d\n", pos+1, var_cnt);
          j=1;
       }
       if(end > pos && !n) {for(i=pos+1; i <= end ; i++) {fprintf(ofp,"%d	ref	0\n", i); fprintf(ofp,"%d	nonref	0\n", i);}}
       if(!j && n) {for(i=beg+1; i <= end ; i++) {fprintf(ofp,"%d	ref	0\n", i); fprintf(ofp,"%d	nonref	0\n", i);}}
       fclose(ofp);
       sprintf(outfn, "%s/Exon%d_depth", pmtr->out_foldr, exn_cnt);
       genecov_draw_plots(pmtr->out_foldr, outfn, exn_cnt);
       fprintf(ofp1,"Exon no = Exon%d\n",exn_cnt);
       fprintf(ofp1,"Genomic coordinates = %s\n",chrpos);
       fprintf(ofp1,"Exon length (bp) = %d\n", esize);
       fprintf(ofp1,"GC%% = %.2f\n", ((float)gc/(float)esize)*100);
       fprintf(ofp1,"Mean depth of coverage = %.2fX\n", ((float)totdepth/(float)esize));
       fprintf(ofp1,"%% covered with at least 1X depth = %.2f\n", ((float)pct1/(float)esize)*100);
       fprintf(ofp1,"%% covered with at least 5X depth = %.2f\n", ((float)pct5/(float)esize)*100);
       fprintf(ofp1,"%% covered with at least 10X depth = %.2f\n", ((float)pct10/(float)esize)*100);
       fprintf(ofp1,"%% covered with at least 20X depth = %.2f\n", ((float)pct20/(float)esize)*100);
       fprintf(ofp1,"%% covered with at least 30X depth = %.2f\n\n", ((float)pct30/(float)esize)*100);
       totaldepth += totdepth;
       //printf("cum depth = %d\n",totaldepth);
       exn_cnt++;
       
}
	
	fprintf(ofp1,"No of exons = %d\n", exn_cnt-1);
	fprintf(ofp1,"Target bed size (bp) = %d\n", totesize);
	fprintf(ofp1,"Mean depth of coverage = %.2fX\n", ((float)totaldepth/(float)(totesize)));
	fclose(ifp);		
	fclose(ofp1);
	
	Arrange_logs(pmtr->out_foldr);
	creathtml(pmtr->out_foldr, exn_cnt-1, totesize, ((float)totaldepth/(float)(totesize)));


	
	
	free(n_plp);
	bam_plp_destroy(mplp);
        bam_index_destroy(idx);
	bam_header_destroy(h);
	bam_close(data->fp);
	bam_iter_destroy(data->iter);
	free(ref);
	free(data);
        return 0;
}

