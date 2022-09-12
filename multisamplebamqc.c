#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bam.h"

void multisamplebamqc_help()
 {
   fprintf(stderr, "\n");
   fprintf(stderr, "Usage: mapinsights multisample-bamqc -i folder_path.list\n\n");
   fprintf(stderr, "About input file structure:\n");
   fprintf(stderr, "folder_path.list file containing a list of bamqc output folder path, one per line.\n\n");
   fprintf(stderr, "There are two tab separated columns in the input file (folder_path.list), \n");
   fprintf(stderr, "first column should be sample_id and the second column should be the path to \n");
   fprintf(stderr, "bamqc output folder. Example structure of the file is given below:\n\n");
   fprintf(stderr, " [sample1	path_to_bamqc_output_folder_of_sample1]\n");
   fprintf(stderr, " [sample2	path_to_bamqc_output_folder_of_sample2]\n");
   fprintf(stderr, " [  .	       .                                      ]\n");
   fprintf(stderr, " [  .	       .                                      ]\n");
   fprintf(stderr, " [sampleN	path_to_bamqc_output_folder_of_sampleN]\n");
   //fprintf(stderr, "	-o output_prefix\n");
   fprintf(stderr, "\n");
   //return 1;
 }

typedef struct {
 char *input_file;
} pramtr;

int main_multisamplebamqc(int argc,char *argv[])
 {
 
      FILE *fp1, *fp2;
   FILE *fp1_t, *fp1_p, *fp2_t, *fp2_p, *fp3_t, *fp3_p;
   FILE *fp4_t, *fp4_p, *fp5_t, *fp5_p, *fp6_t, *fp6_p;
   int i, j, n;
   char ch, sam_id[1024], pth[5000], fl_pth[5100];
   char col1[1000][100], col2[1000][100], outfn[500];
   float column2[12];
   
   pramtr *pmtr = calloc(1, sizeof(pramtr));
   
   while ((n = getopt(argc, argv, "i:h")) >= 0) {
		switch (n) {
			 case 'i': pmtr->input_file = optarg; break; // Input bam file
			 //case 'o': prmtr->out_file = optarg; break;
			 case 'h': multisamplebamqc_help(); return 1; break;
                        }
	}
	if (argc == 1 || optind-1 == argc) {multisamplebamqc_help(); return 1;}
   
   if(access(pmtr->input_file, F_OK) == -1) {printf("Input file does not exist.\n"); return 1;}
   
   fp1 = fopen(pmtr->input_file, "r");
   fp1_t = fopen("GC-content_temp.txt", "w");
   fp1_p = fopen("GC-content_temp-plot.txt", "w");
   fp2_t = fopen("I-Size_temp.txt", "w");
   fp2_p = fopen("I-Size_temp-plot.txt", "w");
   fp3_t = fopen("MapQ_temp.txt", "w");
   fp3_p = fopen("MapQ_temp-plot.txt", "w");
   fp4_t = fopen("Cycle_temp.txt", "w");
   fp4_p = fopen("Cycle_temp-plot.txt", "w");
   fp5_t = fopen("Mismatch-content_temp.txt", "w");
   fp5_p = fopen("Mismatch-content_temp-plot.txt", "w");
   fp6_t = fopen("Substitution-categories_temp.txt", "w");
   fp6_p = fopen("Substitution-categories_temp-plot.txt", "w");
   fprintf(fp1_p,"A	B	samples\n");
   fprintf(fp2_p,"A	B	samples\n");
   fprintf(fp3_p,"A	B	samples\n");
   fprintf(fp4_p,"A	B	samples\n");
   fprintf(fp5_p,"A	B	samples\n");
   fprintf(fp6_p,"A	B	samples\n");
   while((ch=fgetc(fp1)) != EOF)
   {
    //printf("I am in\n");
    i=0;
    sam_id[i]=ch; i++;
    while(((ch=fgetc(fp1)) != '\t'))
    {
      sam_id[i]=ch; i++;
    }
    sam_id[i]='\0';
    //printf("I am in\n");
    i=0;
    while(((ch=fgetc(fp1)) != '\n'))
    {
      pth[i]=ch; i++;
    }
    pth[i]='\0';
    
    printf("%s	%s\n",sam_id, pth);
    //printf("Level1\n");
    sprintf(fl_pth, "%s/logs/Overall_GC_dist.txt", pth);
    //printf("Level2\n");
    if(access(fl_pth, F_OK) == -1) {printf("Unable to locate bamqc output file in given location.\n"); return 1;}
    
    fp2 = fopen(fl_pth, "r");
    //printf("Level3 : %s\n", fl_pth);
    while((ch=fgetc(fp2)) != '\n');
    n=0;
    //printf("Level4\n");
    while((ch=fgetc(fp2)) != EOF){
    //printf("%c",ch);
    i=0;
    col1[n][i]=ch; i++;
    while(((ch=fgetc(fp2)) != '\t'))
    {
      col1[n][i]=ch; i++;
    }
    col1[n][i]='\0';
    i=0;
    
    while(((ch=fgetc(fp2)) != '\n'))
    {
      col2[n][i]=ch; i++;
    }
    col2[n][i]='\0';
    
    //printf("Level5 :  %d	%s	%s\n", n, col1[n], col2[n]);
    n++;
    }
    fclose(fp2);
    //printf("Level6\n");
    fprintf(fp1_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      fprintf(fp1_t, "	%s", col2[i]);
      fprintf(fp1_p, "%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    fprintf(fp1_t,"\n");
    
    
    ////////////Insert-size/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_Insertsize_dist.txt", pth);
    if(access(fl_pth, F_OK) == -1) {printf("Unable to locate bamqc output file in given location.\n"); return 1;}
    
    
    fp2 = fopen(fl_pth, "r");
    //printf("Level3 : %s\n", fl_pth);
    while((ch=fgetc(fp2)) != '\n');
    n=0;
    //printf("Level4\n");
    while((ch=fgetc(fp2)) != EOF){
    //printf("%c",ch);
    i=0;
    col1[n][i]=ch; i++;
    while(((ch=fgetc(fp2)) != '\t'))
    {
      col1[n][i]=ch; i++;
    }
    col1[n][i]='\0';
    i=0;
    
    while(((ch=fgetc(fp2)) != '\n'))
    {
      col2[n][i]=ch; i++;
    }
    col2[n][i]='\0';
    
    //printf("Level5 :  %d	%s	%s\n", n, col1[n], col2[n]);
    n++;
    }
    fclose(fp2);
    //printf("Level6\n");
    fprintf(fp2_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      fprintf(fp2_t, "	%s", col2[i]);
      fprintf(fp2_p, "%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    fprintf(fp2_t,"\n");
    //////////////Insert-size end//////////////////////////
    
    ////////////MapQ/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_Mappingquality_bin_dist.txt", pth);
    if(access(fl_pth, F_OK) == -1) {printf("Unable to locate bamqc output file in given location.\n"); return 1;}
    
    
    fp2 = fopen(fl_pth, "r");
    //printf("Level3 : %s\n", fl_pth);
    while((ch=fgetc(fp2)) != '\n');
    n=0;
    //printf("Level4\n");
    while((ch=fgetc(fp2)) != EOF){
    //printf("%c",ch);
    i=0;
    col1[n][i]=ch; i++;
    while(((ch=fgetc(fp2)) != '\t'))
    {
      col1[n][i]=ch; i++;
    }
    col1[n][i]='\0';
    i=0;
    
    while(((ch=fgetc(fp2)) != '\n'))
    {
      col2[n][i]=ch; i++;
    }
    col2[n][i]='\0';
    
    //printf("Level5 :  %d	%s	%s\n", n, col1[n], col2[n]);
    n++;
    }
    fclose(fp2);
    //printf("Level6\n");
    fprintf(fp3_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      fprintf(fp3_t, "	%s", col2[i]);
      fprintf(fp3_p, "%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    fprintf(fp3_t,"\n");
    //////////////MapQ end//////////////////////////
    
    ////////////Cycle/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_ReadposTotal_SN_MM.txt", pth);
    if(access(fl_pth, F_OK) == -1) {printf("Unable to locate bamqc output file in given location.\n"); return 1;}
    
    
    fp2 = fopen(fl_pth, "r");
    //printf("Level3 : %s\n", fl_pth);
    while((ch=fgetc(fp2)) != '\n');
    n=0;
    //printf("Level4\n");
    while((ch=fgetc(fp2)) != EOF){
    //printf("%c",ch);
    i=0;
    col1[n][i]=ch; i++;
    while(((ch=fgetc(fp2)) != '\t'))
    {
      col1[n][i]=ch; i++;
    }
    col1[n][i]='\0';
    i=0;
    
    while(((ch=fgetc(fp2)) != '\n'))
    {
      col2[n][i]=ch; i++;
    }
    col2[n][i]='\0';
    
    //printf("Level5 :  %d	%s	%s\n", n, col1[n], col2[n]);
    n++;
    }
    fclose(fp2);
    //printf("Level6\n");
    fprintf(fp4_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      fprintf(fp4_t, "	%s", col2[i]);
      fprintf(fp4_p, "p%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    fprintf(fp4_t,"\n");
    //////////////Cycle end//////////////////////////
    
    ////////////Mismatch-content/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_SN_MM.txt", pth);
    if(access(fl_pth, F_OK) == -1) {printf("Unable to locate bamqc output file in given location.\n"); return 1;}
    
    
    fp2 = fopen(fl_pth, "r");
    //printf("Level3 : %s\n", fl_pth);
    while((ch=fgetc(fp2)) != '\n');
    n=0;
    //printf("Level4\n");
    while((ch=fgetc(fp2)) != EOF){
    //printf("%c",ch);
    i=0;
    col1[n][i]=ch; i++;
    while(((ch=fgetc(fp2)) != '\t'))
    {
      col1[n][i]=ch; i++;
    }
    col1[n][i]='\0';
    i=0;
    
    while(((ch=fgetc(fp2)) != '\n'))
    {
      col2[n][i]=ch; i++;
    }
    col2[n][i]='\0';
    
    //printf("Level5- Mismatch-content :  %d	%s	%s\n", n, col1[n], col2[n]);
    n++;
    }
    fclose(fp2);
    //printf("Level6\n");
    fprintf(fp5_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      fprintf(fp5_t, "	%s", col2[i]);
      fprintf(fp5_p, "%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    fprintf(fp5_t,"\n");
    //////////////Mismatch-content end//////////////////////////
    
    ////////////Substitution-categories/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_BaseChange_SN_MM.txt", pth);
    if(access(fl_pth, F_OK) == -1) {printf("Unable to locate bamqc output file in given location.\n"); return 1;}
    
    
    fp2 = fopen(fl_pth, "r");
    //printf("Level3 : %s\n", fl_pth);
    while((ch=fgetc(fp2)) != '\n');
    n=0;
    //printf("Level4\n");
    while((ch=fgetc(fp2)) != EOF){
    //printf("%c",ch);
    i=0;
    col1[n][i]=ch; i++;
    while(((ch=fgetc(fp2)) != '\t'))
    {
      col1[n][i]=ch; i++;
    }
    col1[n][i]='\0';
    i=0;
    
    fscanf(fp2, "%f", &column2[n]);
    ch=fgetc(fp2);
    //printf("Level5 :  %d	%s	%s\n", n, col1[n], col2[n]);
    n++;
    }
    fclose(fp2);
    //printf("Level6\n");
    
    for(i=0; i<n;i++)
    {
            fprintf(fp6_p, "%s	%f	%s\n", col1[i], column2[i], sam_id);
    }
    fprintf(fp6_t,"%s", sam_id);
    if(column2[7] > column2[9]) {fprintf(fp6_t,"	%f", (column2[7] - column2[9]));} else {fprintf(fp6_t,"	%f", (column2[9] - column2[7]));}
    if(column2[2] > column2[4]) {fprintf(fp6_t,"	%f", (column2[2] - column2[4]));} else {fprintf(fp6_t,"	%f", (column2[4] - column2[2]));}
    if(column2[1] > column2[5]) {fprintf(fp6_t,"	%f", (column2[1] - column2[5]));} else {fprintf(fp6_t,"	%f", (column2[5] - column2[1]));}
    if(column2[6] > column2[10]) {fprintf(fp6_t,"	%f", (column2[6] - column2[10]));} else {fprintf(fp6_t,"	%f", (column2[10] - column2[6]));}
    if(column2[0] > column2[3]) {fprintf(fp6_t,"	%f", (column2[0] - column2[3]));} else {fprintf(fp6_t,"	%f", (column2[3] - column2[0]));}
    if(column2[8] > column2[11]) {fprintf(fp6_t,"	%f", (column2[8] - column2[11]));} else {fprintf(fp6_t,"	%f", (column2[11] - column2[8]));}
    fprintf(fp6_t,"\n");
    //////////////Substitution-categories end//////////////////////////
   }
   fclose(fp1);
   
   fclose(fp1_t);
   fclose(fp1_p);
   fclose(fp2_t);
   fclose(fp2_p);
   fclose(fp3_t);
   fclose(fp3_p);
   fclose(fp4_t);
   fclose(fp4_p);
   fclose(fp5_t);
   fclose(fp5_p);
   fclose(fp6_t);
   fclose(fp6_p);
   
   fp1 = fopen("mapinsights_multisample_bamqc_v1.r","w");
   fprintf(fp1,"library(\"ggplot2\")\n");
fprintf(fp1,"library(\"ggdendro\")\n");
fprintf(fp1,"library(\"gridExtra\")\n");
fprintf(fp1,"tbl <- read.table(\"GC-content_temp.txt\", header = FALSE, sep = \"\\t\", row.names=1)\n");
fprintf(fp1,"d <- dist(tbl, method = \"euclidean\")\n");
fprintf(fp1,"hc <- hclust(d, method = \"ward.D\")\n");
fprintf(fp1,"write.table(hc$labels[hc$order], file = \"GC-content_sample-clustering-order.txt\", row.names=FALSE, col.names=FALSE, sep=\"\t\", quote = FALSE)\n");

fprintf(fp1,"dg1 <- ggdendrogram(hc) + ggtitle(\"Clustering based on GC-content distribution\")\n");

fprintf(fp1,"plt <- read.table(\"GC-content_temp-plot.txt\", header = TRUE, sep = \"\t\")\n");

fprintf(fp1,"dg2 <- ggplot(data=plt, aes(x=A, y=B, group=samples)) + geom_line(aes(color = samples)) + xlab(\"GC content (%%)\") + ylab(\"reads(%%)\") + ggtitle(\"All samples - GC-content distribution\")\n");

fprintf(fp1,"gg1 <- ggplot_gtable(ggplot_build(dg1))\n");
fprintf(fp1,"gg2 <- ggplot_gtable(ggplot_build(dg2))\n");
fprintf(fp1,"pdf(\"Mulit-sample_bamqc_output.pdf\")\n");
fprintf(fp1,"grid.arrange(gg2, gg1)\n");




fprintf(fp1,"tbl <- read.table(\"I-Size_temp.txt\", header = FALSE, sep = \"\\t\", row.names=1)\n");
fprintf(fp1,"d <- dist(tbl, method = \"euclidean\")\n");
fprintf(fp1,"hc <- hclust(d, method = \"ward.D\")\n");
fprintf(fp1,"write.table(hc$labels[hc$order], file = \"Insert-size_sample-clustering-order.txt\", row.names=FALSE, col.names=FALSE, sep=\"\t\", quote = FALSE)\n");

fprintf(fp1,"dg1 <- ggdendrogram(hc) + ggtitle(\"Clustering based on insert-size distribution\")\n");

fprintf(fp1,"plt <- read.table(\"I-Size_temp-plot.txt\", header = TRUE, sep = \"\t\")\n");

fprintf(fp1,"dg2 <- ggplot(data=plt, aes(x=A, y=B, group=samples)) + geom_line(aes(color = samples)) + xlab(\"insert size (bp)\") + ylab(\"number of reads\") + ggtitle(\"All samples - insert-size distribution\")\n");

fprintf(fp1,"gg1 <- ggplot_gtable(ggplot_build(dg1))\n");
fprintf(fp1,"gg2 <- ggplot_gtable(ggplot_build(dg2))\n");

fprintf(fp1,"grid.arrange(gg2, gg1)\n");




fprintf(fp1,"tbl <- read.table(\"MapQ_temp.txt\", header = FALSE, sep = \"\\t\", row.names=1)\n");
fprintf(fp1,"d <- dist(tbl, method = \"euclidean\")\n");
fprintf(fp1,"hc <- hclust(d, method = \"ward.D\")\n");
fprintf(fp1,"write.table(hc$labels[hc$order], file = \"Mapping-quality_sample-clustering-order.txt\", row.names=FALSE, col.names=FALSE, sep=\"\t\", quote = FALSE)\n");

fprintf(fp1,"dg1 <- ggdendrogram(hc) + ggtitle(\"Clustering based on mapping-quality distribution\")\n");

fprintf(fp1,"plt <- read.table(\"MapQ_temp-plot.txt\", header = TRUE, sep = \"\t\")\n");
fprintf(fp1,"dg2 <- ggplot(data=plt, aes(x=A, y=B, group=samples)) + geom_bar(stat=\"identity\", position = \"dodge\", aes(color = samples)) + xlab(\"mapping quality bin\") + ylab(\"reads(%%)\") + ggtitle(\"All samples - mapping-quality distribution\")\n");

fprintf(fp1,"gg1 <- ggplot_gtable(ggplot_build(dg1))\n");
fprintf(fp1,"gg2 <- ggplot_gtable(ggplot_build(dg2))\n");

fprintf(fp1,"grid.arrange(gg2, gg1)\n");




fprintf(fp1,"tbl <- read.table(\"Mismatch-content_temp.txt\", header = FALSE, sep = \"\\t\", row.names=1)\n");
fprintf(fp1,"d <- dist(tbl, method = \"euclidean\")\n");
fprintf(fp1,"hc <- hclust(d, method = \"ward.D\")\n");
fprintf(fp1,"write.table(hc$labels[hc$order], file = \"Mismatch-content_sample-clustering-order.txt\", row.names=FALSE, col.names=FALSE, sep=\"\t\", quote = FALSE)\n");

fprintf(fp1,"dg1 <- ggdendrogram(hc) + ggtitle(\"Clustering based on mismatch-content distribution\")\n");

fprintf(fp1,"plt <- read.table(\"Mismatch-content_temp-plot.txt\", header = TRUE, sep = \"\t\")\n");
fprintf(fp1,"mmodr <- c(\"1\",\"2\",\"3\",\"4\",\">=5\")\n");

fprintf(fp1,"dg2 <- ggplot(data=plt, aes(x=A, y=B, group=samples)) + geom_bar(stat=\"identity\", position = \"dodge\", aes(color = samples)) + scale_x_discrete(limits = mmodr) + xlab(\"mismatch counts\") + ylab(\"reads(%%)\") + ggtitle(\"All samples - mismatch-content distribution\")\n");
fprintf(fp1,"gg1 <- ggplot_gtable(ggplot_build(dg1))\n");
fprintf(fp1,"gg2 <- ggplot_gtable(ggplot_build(dg2))\n");

fprintf(fp1,"grid.arrange(gg2, gg1)\n");



fprintf(fp1,"tbl <- read.table(\"Substitution-categories_temp.txt\", header = FALSE, sep = \"\\t\", row.names=1)\n");
fprintf(fp1,"d <- dist(tbl, method = \"euclidean\")\n");
fprintf(fp1,"hc <- hclust(d, method = \"ward.D\")\n");
fprintf(fp1,"write.table(hc$labels[hc$order], file = \"Substitution-categories_sample-clustering-order.txt\", row.names=FALSE, col.names=FALSE, sep=\"\t\", quote = FALSE)\n");

fprintf(fp1,"dg1 <- ggdendrogram(hc) + ggtitle(\"Clustering based on distribution of substitution-categories\")\n");

fprintf(fp1,"plt <- read.table(\"Substitution-categories_temp-plot.txt\", header = TRUE, sep = \"\t\")\n");
fprintf(fp1,"morder <- c(\"C>T\",\"G>A\",\"A>G\",\"T>C\",\"A>C\",\"T>G\",\"C>A\",\"G>T\",\"A>T\",\"T>A\",\"C>G\",\"G>C\")\n");

fprintf(fp1,"dg2 <- ggplot(plt, aes(x=A, y=B)) + geom_boxplot(fill='orange', color=\"gray30\") + geom_jitter(position=position_jitter(0.2)) + scale_x_discrete(limits = morder) + xlab(\"change\") + ylab(\"change(%%)\") + ggtitle(\"All samples - distribution of substitution-categories\")\n");

fprintf(fp1,"gg1 <- ggplot_gtable(ggplot_build(dg1))\n");
fprintf(fp1,"gg2 <- ggplot_gtable(ggplot_build(dg2))\n");

fprintf(fp1,"grid.arrange(gg2, gg1)\n");



fprintf(fp1,"tbl <- read.table(\"Cycle_temp.txt\", header = FALSE, sep = \"\\t\", row.names=1, col.names = paste0(\"V\",seq_len(max(count.fields(\"Cycle_temp.txt\", sep = \"\t\")))), fill = TRUE)\n");
fprintf(fp1,"d <- dist(tbl, method = \"euclidean\")\n");
fprintf(fp1,"hc <- hclust(d, method = \"ward.D\")\n");
fprintf(fp1,"write.table(hc$labels[hc$order], file = \"Per-cycle_sample-clustering-order.txt\", row.names=FALSE, col.names=FALSE, sep=\"\t\", quote = FALSE)\n");

fprintf(fp1,"dg1 <- ggdendrogram(hc) + ggtitle(\"Clustering based on per-cycle substitution rate\")\n");

fprintf(fp1,"plt <- read.table(\"Cycle_temp-plot.txt\", header = TRUE, sep = \"\t\")\n");
fprintf(fp1,"plt$A <- factor(plt$A, levels = unique(plt$A))\n");

fprintf(fp1,"dg2 <- ggplot(plt, aes(x=A, y=B)) + geom_boxplot(fill='deepskyblue1') + theme(axis.text.x = element_text(angle = 90)) + xlab(\"read position (bp)\") + ylab(\"reads having mismatches (%%)\") + ggtitle(\"All samples - distribution of per-cycle substitution rate\")\n");
fprintf(fp1,"gg1 <- ggplot_gtable(ggplot_build(dg1))\n");
fprintf(fp1,"gg2 <- ggplot_gtable(ggplot_build(dg2))\n");

fprintf(fp1,"grid.arrange(gg2, gg1)\n");

fprintf(fp1,"dev.off()\n");

fclose(fp1);
      int systemRet = system("R CMD BATCH mapinsights_multisample_bamqc_v1.r");
if(systemRet == -1){printf("Fail\n\n");}
      systemRet = system("rm -f mapinsights_multisample_bamqc_v1.* *.RData");
if(systemRet == -1){printf("Fail\n\n");}
      systemRet = system("rm -f Rplots.pdf Cycle_temp.txt I-Size_temp.txt Mismatch-content_temp.txt GC-content_temp.txt MapQ_temp.txt Substitution-categories_temp.txt");
if(systemRet == -1){printf("Fail\n\n");}
      systemRet = system("rm -f Cycle_temp-plot.txt I-Size_temp-plot.txt Mismatch-content_temp-plot.txt GC-content_temp-plot.txt MapQ_temp-plot.txt Substitution-categories_temp-plot.txt");
if(systemRet == -1){printf("Fail\n\n");}
      
       
}

