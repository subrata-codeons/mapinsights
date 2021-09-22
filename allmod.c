#include <stdio.h>
#include <string.h>

char *MAPINSIGHTS_VERSION = "1.0";

int main_bamqc(int argc,char *argv[]);
int main_genedepth(int argc, char *argv[]);
int main_siteinfo(int argc, char *argv[]);
int main_jumpreads(int argc, char *argv[]);

int mapinsights_help()
 {
   fprintf(stderr, "\n");
   fprintf(stderr, "Program: mapinsights\n");
   fprintf(stderr, "Version: %s\n\n", MAPINSIGHTS_VERSION);
   fprintf(stderr, "Usage:   mapinsights <command> [options]\n\n");
   fprintf(stderr, "Command: bamqc		QC of alignment file\n");
   fprintf(stderr, "         genedepth	Estimate exon-wise bed coverage\n");
   fprintf(stderr, "         siteinfo	Details information about genomic site(s)\n");
   fprintf(stderr, "         jumpreads	Extract reads with jump alignment\n");
   fprintf(stderr, "\n");
   return 1;
   }

int main(int argc, char *argv[])
{
  //int n=0;
  //while ((n = getopt(argc, argv, "h")) >= 0) {switch (n) {case 'h': mapinsights_help(); break;}}
  //printf("passed command : %s\n",argv[1]);	
  if (argc < 2) return mapinsights_help();
  if (strcmp(argv[1], "bamqc") == 0) return main_bamqc(argc-1, argv+1);
  else if (strcmp(argv[1], "genedepth") == 0) return main_genedepth(argc-1, argv+1);
  else if (strcmp(argv[1], "siteinfo") == 0) return main_siteinfo(argc-1, argv+1);
  else if (strcmp(argv[1], "jumpreads") == 0) return main_jumpreads(argc-1, argv+1);
  //else if (n != 0) {printf("Unknown command\n"); mapinsights_help();}

  return 0;
}
