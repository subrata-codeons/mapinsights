#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bam.h"

void batchplotbamqc_help()
 {
   fprintf(stderr, "\n");
   fprintf(stderr, "Usage: mapinsights multisample-bamqc -1 folder_path_for_batch1.list -2 folder_path_for_batch2.list\n\n");
   fprintf(stderr, "About input file structure:\n");
   fprintf(stderr, "folder_path_for_batch1/2.list file containing a list of bamqc output folder path, one per line.\n\n");
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
 char *input_file1, *input_file2;
} pramtr;

int main_batchplotbamqc(int argc,char *argv[])
{
   FILE *fp1, *fp2;
   FILE *fp1_b1, *fp1_b2, *fp2_b1, *fp2_b2, *fp3_b1, *fp3_b2;
   FILE *fp4_b1, *fp4_b2, *fp5_b1, *fp5_b2, *fp6_b1, *fp6_b2;
   int i, j, n;
   char ch, sam_id[1024], pth[5000], fl_pth[5100];
   char col1[1000][100], col2[1000][100], outfn[500];
   float column2[12];
   
   pramtr *pmtr = calloc(1, sizeof(pramtr));
  
   while ((n = getopt(argc, argv, "1:2:h")) >= 0) {
		switch (n) {
			 case '1': pmtr->input_file1 = optarg; break; 
			 case '2': pmtr->input_file2 = optarg; break;
			 case 'h': batchplotbamqc_help(); return 1; break;
                        }
	}
	if (argc == 1 || optind-1 == argc) {batchplotbamqc_help(); return 1;}
   
   if(access(pmtr->input_file1, F_OK) == -1) {printf("Input file does not exist.\n"); return 1;}
   if(access(pmtr->input_file2, F_OK) == -1) {printf("Input file does not exist.\n"); return 1;}
   
   fp1 = fopen(pmtr->input_file1, "r");
   fp1_b1 = fopen("GC-content_temp-plot_b1.txt", "w");
   fp1_b2 = fopen("GC-content_temp-plot_b2.txt", "w");
   fp2_b1 = fopen("I-Size_temp-plot_b1.txt", "w");
   fp2_b2 = fopen("I-Size_temp-plot_b2.txt", "w");
   fp3_b1 = fopen("MapQ_temp-plot_b1.txt", "w");
   fp3_b2 = fopen("MapQ_temp-plot_b2.txt", "w");
   fp4_b1 = fopen("Cycle_temp-plot_b1.txt", "w");
   fp4_b2 = fopen("Cycle_temp-plot_b2.txt", "w");
   fp5_b1 = fopen("Mismatch-content_temp-plot_b1.txt", "w");
   fp5_b2 = fopen("Mismatch-content_temp-plot_b2.txt", "w");
   fp6_b1 = fopen("Substitution-categories_temp-plot_b1.txt", "w");
   fp6_b2 = fopen("Substitution-categories_temp-plot_b2.txt", "w");
   fprintf(fp1_b1,"A	B	samples\n");
   fprintf(fp2_b1,"A	B	samples\n");
   fprintf(fp3_b1,"A	B	samples\n");
   fprintf(fp4_b1,"A	B	samples\n");
   fprintf(fp5_b1,"A	B	samples\n");
   fprintf(fp6_b1,"A	B	samples\n");
   
   fprintf(fp1_b2,"A	B	samples\n");
   fprintf(fp2_b2,"A	B	samples\n");
   fprintf(fp3_b2,"A	B	samples\n");
   fprintf(fp4_b2,"A	B	samples\n");
   fprintf(fp5_b2,"A	B	samples\n");
   fprintf(fp6_b2,"A	B	samples\n");
   
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
    //fprintf(fp1_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      //fprintf(fp1_t, "	%s", col2[i]);
      fprintf(fp1_b1, "%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    //fprintf(fp1_t,"\n");
    
    
    ////////////Insert-size/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_Insertsize_dist.txt", pth);
    //printf("Level2\n");
    
    
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
    //fprintf(fp2_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      //fprintf(fp2_t, "	%s", col2[i]);
      fprintf(fp2_b1, "%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    //fprintf(fp2_t,"\n");
    //////////////Insert-size end//////////////////////////
    
    ////////////MapQ/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_Mappingquality_bin_dist.txt", pth);
    //printf("Level2\n");
    
    
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
    //fprintf(fp3_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      //fprintf(fp3_t, "	%s", col2[i]);
      fprintf(fp3_b1, "%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    //fprintf(fp3_t,"\n");
    //////////////MapQ end//////////////////////////
    
    ////////////Cycle/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_ReadposTotal_SN_MM.txt", pth);
    //printf("Level2\n");
    
    
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
    //fprintf(fp4_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      //fprintf(fp4_t, "	%s", col2[i]);
      fprintf(fp4_b1, "p%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    //fprintf(fp4_t,"\n");
    //////////////Cycle end//////////////////////////
    
    ////////////Mismatch-content/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_SN_MM.txt", pth);
    //printf("Level2\n");
    
    
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
    //fprintf(fp5_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      //fprintf(fp5_t, "	%s", col2[i]);
      fprintf(fp5_b1, "%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    //fprintf(fp5_t,"\n");
    //////////////Mismatch-content end//////////////////////////
    
    ////////////Substitution-categories/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_BaseChange_SN_MM.txt", pth);
    //printf("Level2\n");
    
    
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
            fprintf(fp6_b1, "%s	%f	%s\n", col1[i], column2[i], sam_id);
    }
    //fprintf(fp6_t,"%s", sam_id);
    //if(column2[7] > column2[9]) {fprintf(fp6_t,"	%f", (column2[7] - column2[9]));} else {fprintf(fp6_t,"	%f", (column2[9] - column2[7]));}
    //if(column2[2] > column2[4]) {fprintf(fp6_t,"	%f", (column2[2] - column2[4]));} else {fprintf(fp6_t,"	%f", (column2[4] - column2[2]));}
    //if(column2[1] > column2[5]) {fprintf(fp6_t,"	%f", (column2[1] - column2[5]));} else {fprintf(fp6_t,"	%f", (column2[5] - column2[1]));}
    //if(column2[6] > column2[10]) {fprintf(fp6_t,"	%f", (column2[6] - column2[10]));} else {fprintf(fp6_t,"	%f", (column2[10] - column2[6]));}
    //if(column2[0] > column2[3]) {fprintf(fp6_t,"	%f", (column2[0] - column2[3]));} else {fprintf(fp6_t,"	%f", (column2[3] - column2[0]));}
    //if(column2[8] > column2[11]) {fprintf(fp6_t,"	%f", (column2[8] - column2[11]));} else {fprintf(fp6_t,"	%f", (column2[11] - column2[8]));}
    //fprintf(fp6_t,"\n");
    //////////////Substitution-categories end//////////////////////////
   }
   fclose(fp1);
   
   fclose(fp1_b1);
   fclose(fp2_b1);
   fclose(fp3_b1);
   fclose(fp4_b1);
   fclose(fp5_b1);
   fclose(fp6_b1);
      
   /////////////Batch2 start//////////////////
   fp1 = fopen(pmtr->input_file2, "r");
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
    //fprintf(fp1_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      //fprintf(fp1_t, "	%s", col2[i]);
      fprintf(fp1_b2, "%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    //fprintf(fp1_t,"\n");
    
    
    ////////////Insert-size/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_Insertsize_dist.txt", pth);
    //printf("Level2\n");
    
    
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
    //fprintf(fp2_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      //fprintf(fp2_t, "	%s", col2[i]);
      fprintf(fp2_b2, "%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    //fprintf(fp2_t,"\n");
    //////////////Insert-size end//////////////////////////
    
    ////////////MapQ/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_Mappingquality_bin_dist.txt", pth);
    //printf("Level2\n");
    
    
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
    //fprintf(fp3_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      //fprintf(fp3_t, "	%s", col2[i]);
      fprintf(fp3_b2, "%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    //fprintf(fp3_t,"\n");
    //////////////MapQ end//////////////////////////
    
    ////////////Cycle/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_ReadposTotal_SN_MM.txt", pth);
    //printf("Level2\n");
    
    
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
    //fprintf(fp4_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      //fprintf(fp4_t, "	%s", col2[i]);
      fprintf(fp4_b2, "p%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    //fprintf(fp4_t,"\n");
    //////////////Cycle end//////////////////////////
    
    ////////////Mismatch-content/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_SN_MM.txt", pth);
    //printf("Level2\n");
    
    
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
    //fprintf(fp5_t,"%s", sam_id);
    for(i=0; i<n;i++)
    {
      //fprintf(fp5_t, "	%s", col2[i]);
      fprintf(fp5_b2, "%s	%s	%s\n", col1[i], col2[i], sam_id);
    }
    //fprintf(fp5_t,"\n");
    //////////////Mismatch-content end//////////////////////////
    
    ////////////Substitution-categories/////////////////
    fl_pth[0]='\0';
    sprintf(fl_pth, "%s/logs/Overall_BaseChange_SN_MM.txt", pth);
    //printf("Level2\n");
    
    
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
            fprintf(fp6_b2, "%s	%f	%s\n", col1[i], column2[i], sam_id);
    }
    //fprintf(fp6_t,"%s", sam_id);
    //if(column2[7] > column2[9]) {fprintf(fp6_t,"	%f", (column2[7] - column2[9]));} else {fprintf(fp6_t,"	%f", (column2[9] - column2[7]));}
    //if(column2[2] > column2[4]) {fprintf(fp6_t,"	%f", (column2[2] - column2[4]));} else {fprintf(fp6_t,"	%f", (column2[4] - column2[2]));}
    //if(column2[1] > column2[5]) {fprintf(fp6_t,"	%f", (column2[1] - column2[5]));} else {fprintf(fp6_t,"	%f", (column2[5] - column2[1]));}
    //if(column2[6] > column2[10]) {fprintf(fp6_t,"	%f", (column2[6] - column2[10]));} else {fprintf(fp6_t,"	%f", (column2[10] - column2[6]));}
    //if(column2[0] > column2[3]) {fprintf(fp6_t,"	%f", (column2[0] - column2[3]));} else {fprintf(fp6_t,"	%f", (column2[3] - column2[0]));}
    //if(column2[8] > column2[11]) {fprintf(fp6_t,"	%f", (column2[8] - column2[11]));} else {fprintf(fp6_t,"	%f", (column2[11] - column2[8]));}
    //fprintf(fp6_t,"\n");
    //////////////Substitution-categories end//////////////////////////
   }
   fclose(fp1);
   
   fclose(fp1_b2);
   fclose(fp2_b2);
   fclose(fp3_b2);
   fclose(fp4_b2);
   fclose(fp5_b2);
   fclose(fp6_b2);
   /////////////Batch2 end////////////////////
   fp1 = fopen("mapinsights_multisample_batchplot_v1.r","w");
fprintf(fp1, "library(\"ggplot2\")\n");
fprintf(fp1, "library(\"gridExtra\")\n");


fprintf(fp1, "plt <- read.table(\"GC-content_temp-plot_b1.txt\", header = TRUE, sep = \"\t\")\n");

fprintf(fp1, "dg1 <- ggplot(data=plt, aes(x=A, y=B, group=samples)) + geom_line(aes(color = samples)) + xlab(\"GC content (%%)\") + ylab(\"pct of reads\") + ggtitle(\"All samples (batch1) - GC-content distribution\")\n");

fprintf(fp1, "plt <- read.table(\"GC-content_temp-plot_b2.txt\", header = TRUE, sep = \"\t\")\n");

fprintf(fp1, "dg2 <- ggplot(data=plt, aes(x=A, y=B, group=samples)) + geom_line(aes(color = samples)) + xlab(\"GC content (%%)\") + ylab(\"pct of reads\") + ggtitle(\"All samples (batch2) - GC-content distribution\")\n");

fprintf(fp1, "gg1 <- ggplot_gtable(ggplot_build(dg1))\n");
fprintf(fp1, "gg2 <- ggplot_gtable(ggplot_build(dg2))\n");
fprintf(fp1, "pdf(\"Mulit-sample_batch-plots.pdf\")\n");
fprintf(fp1, "grid.arrange(gg1, gg2)\n");

//#fprintf(fp1, "dev.off()\n");


fprintf(fp1, "plt <- read.table(\"I-Size_temp-plot_b1.txt\", header = TRUE, sep = \"\t\")\n");

fprintf(fp1, "dg1 <- ggplot(data=plt, aes(x=A, y=B, group=samples)) + geom_line(aes(color = samples)) + xlab(\"insert size (bp)\") + ylab(\"no of reads\") + ggtitle(\"All samples (batch1) - insert-size distribution\")\n");

fprintf(fp1, "plt <- read.table(\"I-Size_temp-plot_b2.txt\", header = TRUE, sep = \"\t\")\n");

fprintf(fp1, "dg2 <- ggplot(data=plt, aes(x=A, y=B, group=samples)) + geom_line(aes(color = samples)) + xlab(\"insert size (bp)\") + ylab(\"no of reads\") + ggtitle(\"All samples (batch2) - insert-size distribution\")\n");

fprintf(fp1, "gg1 <- ggplot_gtable(ggplot_build(dg1))\n");
fprintf(fp1, "gg2 <- ggplot_gtable(ggplot_build(dg2))\n");
//#fprintf(fp1, "pdf(\"plots.pdf\")\n");
fprintf(fp1, "grid.arrange(gg1, gg2)\n");

//#fprintf(fp1, "dev.off()\n");




fprintf(fp1, "plt <- read.table(\"MapQ_temp-plot_b1.txt\", header = TRUE, sep = \"\t\")\n");
fprintf(fp1, "dg1 <- ggplot(data=plt, aes(x=A, y=B, group=samples)) + geom_bar(stat=\"identity\", position = \"dodge\", aes(color = samples)) + xlab(\"mapping quality bin\") + ylab(\"pct of reads\") + ggtitle(\"All samples (batch1) - mapping-quality distribution\")\n");

fprintf(fp1, "plt <- read.table(\"MapQ_temp-plot_b2.txt\", header = TRUE, sep = \"\t\")\n");
fprintf(fp1, "dg2 <- ggplot(data=plt, aes(x=A, y=B, group=samples)) + geom_bar(stat=\"identity\", position = \"dodge\", aes(color = samples)) + xlab(\"mapping quality bin\") + ylab(\"pct of reads\") + ggtitle(\"All samples (batch2) - mapping-quality distribution\")\n");

fprintf(fp1, "gg1 <- ggplot_gtable(ggplot_build(dg1))\n");
fprintf(fp1, "gg2 <- ggplot_gtable(ggplot_build(dg2))\n");
//#fprintf(fp1, "pdf(\"plots.pdf\")\n");
fprintf(fp1, "grid.arrange(gg1, gg2)\n");

//#fprintf(fp1, "dev.off()\n");




fprintf(fp1, "plt <- read.table(\"Mismatch-content_temp-plot_b1.txt\", header = TRUE, sep = \"\t\")\n");
fprintf(fp1, "mmodr <- c(\"1\",\"2\",\"3\",\"4\",\">=5\")\n");

fprintf(fp1, "dg1 <- ggplot(data=plt, aes(x=A, y=B, group=samples)) + geom_bar(stat=\"identity\", position = \"dodge\", aes(color = samples)) + scale_x_discrete(limits = mmodr) + xlab(\"mismatch counts\") + ylab(\"pct of reads\") + ggtitle(\"All samples (batch1) - mismatch-content distribution\")\n");

fprintf(fp1, "plt <- read.table(\"Mismatch-content_temp-plot_b2.txt\", header = TRUE, sep = \"\t\")\n");

fprintf(fp1, "dg2 <- ggplot(data=plt, aes(x=A, y=B, group=samples)) + geom_bar(stat=\"identity\", position = \"dodge\", aes(color = samples)) + scale_x_discrete(limits = mmodr) + xlab(\"mismatch counts\") + ylab(\"pct of reads\") + ggtitle(\"All samples (batch2) - mismatch-content distribution\")\n");

fprintf(fp1, "gg1 <- ggplot_gtable(ggplot_build(dg1))\n");
fprintf(fp1, "gg2 <- ggplot_gtable(ggplot_build(dg2))\n");
//#fprintf(fp1, "pdf(\"plots.pdf\")\n");
fprintf(fp1, "grid.arrange(gg1, gg2)\n");

//#fprintf(fp1, "dev.off()\n");




fprintf(fp1, "plt <- read.table(\"Substitution-categories_temp-plot_b1.txt\", header = TRUE, sep = \"\t\")\n");
fprintf(fp1, "morder <- c(\"C>T\",\"G>A\",\"A>G\",\"T>C\",\"A>C\",\"T>G\",\"C>A\",\"G>T\",\"A>T\",\"T>A\",\"C>G\",\"G>C\")\n");

fprintf(fp1, "dg1 <- ggplot(plt, aes(x=A, y=B)) + geom_boxplot(fill='orange', color=\"gray30\") + geom_jitter(position=position_jitter(0.2)) + scale_x_discrete(limits = morder) + xlab(\"change\") + ylab(\"pct of change\") + ggtitle(\"All samples (batch1) - distribution of substitution-categories\")\n");

fprintf(fp1, "plt <- read.table(\"Substitution-categories_temp-plot_b2.txt\", header = TRUE, sep = \"\t\")\n");

fprintf(fp1, "dg2 <- ggplot(plt, aes(x=A, y=B)) + geom_boxplot(fill='orange', color=\"gray30\") + geom_jitter(position=position_jitter(0.2)) + scale_x_discrete(limits = morder) + xlab(\"change\") + ylab(\"pct of change\") + ggtitle(\"All samples (batch2) - distribution of substitution-categories\")\n");

fprintf(fp1, "gg1 <- ggplot_gtable(ggplot_build(dg1))\n");
fprintf(fp1, "gg2 <- ggplot_gtable(ggplot_build(dg2))\n");
//#fprintf(fp1, "pdf(\"plots.pdf\")\n");
fprintf(fp1, "grid.arrange(gg1, gg2)\n");

//#fprintf(fp1, "dev.off()\n");


fprintf(fp1, "plt <- read.table(\"Cycle_temp-plot_b1.txt\", header = TRUE, sep = \"\t\")\n");
fprintf(fp1, "plt$A <- factor(plt$A, levels = unique(plt$A))\n");

fprintf(fp1, "dg1 <- ggplot(plt, aes(x=A, y=B)) + geom_boxplot(fill='deepskyblue1') + theme(axis.text.x = element_text(angle = 90)) + xlab(\"read position (bp)\") + ylab(\"reads having mismatches (%%)\") + ggtitle(\"All samples (batch1) - distribution of per-cycle substitution rate\")\n");

fprintf(fp1, "plt <- read.table(\"Cycle_temp-plot_b2.txt\", header = TRUE, sep = \"\t\")\n");
fprintf(fp1, "plt$A <- factor(plt$A, levels = unique(plt$A))\n");

fprintf(fp1, "dg2 <- ggplot(plt, aes(x=A, y=B)) + geom_boxplot(fill='deepskyblue1') + theme(axis.text.x = element_text(angle = 90)) + xlab(\"read position (bp)\") + ylab(\"reads having mismatches (%%)\") + ggtitle(\"All samples (batch2) - distribution of per-cycle substitution rate\")\n");

fprintf(fp1, "gg1 <- ggplot_gtable(ggplot_build(dg1))\n");
fprintf(fp1, "gg2 <- ggplot_gtable(ggplot_build(dg2))\n");
//#fprintf(fp1, "pdf(\"plots.pdf\")\n");
fprintf(fp1, "grid.arrange(gg1, gg2)\n");

fprintf(fp1, "dev.off()\n");

fclose(fp1);
      int systemRet = system("R CMD BATCH mapinsights_multisample_batchplot_v1.r");
if(systemRet == -1){printf("Fail\n\n");}
      systemRet = system("rm -f mapinsights_multisample_batchplot_v1.* *.RData");
if(systemRet == -1){printf("Fail\n\n");}
      systemRet = system("rm -f Rplots.pdf GC-content_temp-plot_b1.txt GC-content_temp-plot_b2.txt I-Size_temp-plot_b1.txt I-Size_temp-plot_b2.txt MapQ_temp-plot_b1.txt MapQ_temp-plot_b2.txt Cycle_temp-plot_b1.txt");
if(systemRet == -1){printf("Fail\n\n");}
      systemRet = system("rm -f Cycle_temp-plot_b2.txt Mismatch-content_temp-plot_b1.txt Mismatch-content_temp-plot_b2.txt Substitution-categories_temp-plot_b1.txt Substitution-categories_temp-plot_b2.txt");
if(systemRet == -1){printf("Fail\n\n");}
      
}

