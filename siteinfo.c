
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
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

 void siteinfo_help()
 {
   fprintf(stderr, "\n");
   fprintf(stderr, "Usage: mapinsights siteinfo -r <ref.fa> -o <output-file> -i <aligned.bam>\n\n");
   fprintf(stderr, "Options:\n");
   fprintf(stderr, "	-r  ref.fa      reference fasta\n");
   fprintf(stderr, "	-p  query positions file, one per line [example :chr#	pos]\n");
   fprintf(stderr, "	-i  input file   Alignment file (BAM)\n");
   fprintf(stderr, "	-o  output file name \n");
   fprintf(stderr, "	-s  position        [chr#:position-position]\n");
   fprintf(stderr, "	-h  help\n");
   fprintf(stderr, "\n");
   //return 1;
 }
 void prntheader()
 {
   
   printf("# For each input site, there are two blocks in the output. The first output block provides information like genomic coordinates, reference-base, 20bp flanks with GC%% and depth. Depending on the depth the second block provides following information as listed\n");
   printf("# 1st column  : Aligned bases at input position\n");
   printf("# 2nd column  : Quality value of aligned bases at input position\n");
   printf("# 3rd column  : Strand \n");
   printf("# 4th column  : First in pair [rd1] / Second in pair [rd2]\n");
   printf("# 5th column  : Position in read\n");
   printf("# 6th column  : Proper pair [pp] / Not proper pair [npp]\n");
   printf("# 7th column  : RG-tag @RG:ID\n");
   printf("# 8th column  : MD-tag\n");
   printf("# 9th column  : Number of CIGAR operation\n");
   printf("# 10th column : Mapping quality\n");
   printf("# 11th column : Insertsize\n");
   printf("# 12th column : Readid\n");
 }
typedef struct {
 char *ref_file, *pos_file, *bam_file;
 char *out_file, *site;
} pramtr;
int main_siteinfo(int argc, char *argv[])
{
	int i, n,j, tid, beg, end, pos, *n_plp;
	const bam_pileup1_t *q,*p;
	bam_header_t *h = 0; 
	aux_t *data;
        faidx_t *fai;
	bam_plp_t mplp;
        int tid0 = -1, ref_tid = -1, ref_len;
        char *ref, ch, tchr[20],tpos[50], chrpos[200];
        FILE *ifp;
	FILE *ofp;
	pramtr *pmtr = calloc(1, sizeof(pramtr));
        mplp = calloc(1, sizeof(bam_plp_t));
        while ((n = getopt(argc, argv, "r:p:i:o:s:h")) >= 0) {
             switch (n) {
                    case 'r': pmtr->ref_file = optarg; break; // reference file
                    case 'p': pmtr->pos_file = optarg; break; // bed file
                    case 'i': pmtr->bam_file = optarg; break; // Input bam file
                    case 'o': pmtr->out_file = optarg; break; // output folder path
                    case 's': pmtr->site = optarg; break;
                    case 'h': siteinfo_help(); return 1; break;
                    }
        }
        if (argc == 1 || optind-1 == argc) {siteinfo_help(); return 1;}

        if(access(pmtr->ref_file, F_OK) == -1) {printf("Reference file does not exist.\n"); return 1;}

        if(access(pmtr->bam_file, F_OK) == -1) {printf("Alignment(.bam) file does not exist.\n"); return 1;}

        if(pmtr->pos_file) {if(access(pmtr->pos_file, F_OK) == -1) {printf("Position file does not exist.\n"); return 1;}}
	
	if(pmtr->pos_file && pmtr->site) {printf("Please provide either -s or -p option, not both.\n"); return 1;}
	
        if(pmtr->out_file) 
        {
           ofp = fopen(pmtr->out_file, "w");
           if(ofp == NULL) {printf("Unable to open output file\n"); return 1;}
        }
	
	if(pmtr->out_file) {
	fprintf(ofp,"# mapinsights siteinfo output\n");
	fprintf(ofp,"# mapinsights version-1.0.0\n");
	fprintf(ofp,"#Command = mapinsights");
	for(i=0;i<argc; i++) {fprintf(ofp," %s", argv[i]);}
	fprintf(ofp,"\n\n#Alignment file = %s\n", pmtr->bam_file);
	fprintf(ofp,"#Reference file = %s\n", pmtr->ref_file);
	fprintf(ofp,"#Output file = %s\n", pmtr->out_file);
	if(pmtr->pos_file) {fprintf(ofp,"#Position file = %s\n\n", pmtr->pos_file);}
	else if(pmtr->site) {fprintf(ofp,"#Coordinate = %s\n\n", pmtr->site);}
	fprintf(ofp,"# For each input site, there are two blocks in the output. The first output block provides information like genomic coordinates, reference-base, 20bp flanks with GC%% and depth. Depending on the depth the second block provides following information as listed\n");
   fprintf(ofp,"# 1st column  : Aligned bases at input position\n");
   fprintf(ofp,"# 2nd column  : Quality value of aligned bases at input position\n");
   fprintf(ofp,"# 3rd column  : Strand \n");
   fprintf(ofp,"# 4th column  : First in pair [rd1] / Second in pair [rd2]\n");
   fprintf(ofp,"# 5th column  : Position in read\n");
   fprintf(ofp,"# 6th column  : Proper pair [pp] / Not proper pair [npp]\n");
   fprintf(ofp,"# 7th column  : RG-tag @RG:ID\n");
   fprintf(ofp,"# 8th column  : MD-tag\n");
   fprintf(ofp,"# 9th column  : Number of CIGAR operation\n");
   fprintf(ofp,"# 10th column : Mapping quality\n");
   fprintf(ofp,"# 11th column : Insertsize\n");
   fprintf(ofp,"# 12th column : Readid\n");
	
	}
	else
	{
	printf("# mapinsights siteinfo output\n");
	printf("# mapinsights version-1.0.0\n");
	printf("#Command = mapinsights");
	for(i=0;i<argc; i++) {printf(" %s", argv[i]);}
	printf("\n\n#Alignment file = %s\n", pmtr->bam_file);
	printf("#Reference file = %s\n", pmtr->ref_file);
	if(pmtr->pos_file) {printf("#Position file = %s\n\n", pmtr->pos_file);}
	else if(pmtr->site) {printf("#Coordinate = %s\n\n", pmtr->site);}
        prntheader();
        }
        //n = argc - optind; // the number of BAMs on the command line
	//n=1;
	data = calloc(1, sizeof(void*)); // data[i] for the i-th input
	beg = 0; end = 1<<30; tid = -1;  // set the default region
	        //chrpos=calloc(500,1);
		bam_header_t *htmp;
               
		data = calloc(1, sizeof(aux_t));
		
		data->fp = bam_open(pmtr->bam_file, "r");
		htmp = bam_header_read(data->fp);         
                fai = fai_load(pmtr->ref_file);
                tid0 = tid;
        if (tid0 >= 0 && fai) { 
                ref = faidx_fetch_seq(fai, htmp->target_name[tid0], 0, 0x7fffffff, &ref_len);
                ref_tid = tid0;

        } else ref_tid = -1, ref = 0;
                bam_index_t *idx = bam_index_load(pmtr->bam_file);
		h = htmp; 
                if(pmtr->site)
                {
                  bam_parse_region(h, pmtr->site, &tid, &beg, &end);


		if (tid >= 0) 
		{ 
			
			data->iter = bam_iter_query(idx, tid, beg, end); 
			 
                        
		
	
        if (tid != ref_tid) 
                        free(ref); ref = 0;
                       if (fai) ref = faidx_fetch_seq(fai,  htmp->target_name[tid], 0, 0x7fffffff, &ref_len);
                        ref_tid = tid;

                }

        if(pmtr->out_file) { fprintf(ofp, "%s    ref_base = %c	",pmtr->site, toupper(ref[beg]));}
        else {printf("%s    ref_base = %c	",pmtr->site, toupper(ref[beg]));}
        n=0;
        if(pmtr->out_file) {
        for(i=10;i>0;i--) {fprintf(ofp,"%c",toupper(ref[beg-i])); if(toupper(ref[beg-i])=='C' || toupper(ref[beg-i])=='G') {n++;}}
        fprintf(ofp," %c ",toupper(ref[beg]));
        for(i=1;i<=10;i++) {fprintf(ofp,"%c",toupper(ref[beg+i])); if(toupper(ref[beg+i])=='C' || toupper(ref[beg+i])=='G') {n++;}}
        fprintf(ofp,"	GC%%[in 20bp flanks] = %.2f	", ((float)n/20.0)*100);
        }
        else
        {
        for(i=10;i>0;i--) {printf("%c",toupper(ref[beg-i])); if(toupper(ref[beg-i])=='C' || toupper(ref[beg-i])=='G') {n++;}}
        printf(" %c ",toupper(ref[beg]));
        for(i=1;i<=10;i++) {printf("%c",toupper(ref[beg+i])); if(toupper(ref[beg+i])=='C' || toupper(ref[beg+i])=='G') {n++;}}
        printf("	GC%%[in 20bp flanks] = %.2f	", ((float)n/20.0)*100);
        }
	mplp = bam_plp_init(read_bam, (void*)data); 
	n_plp = calloc(1, sizeof(int)); 

        j=0;
	while((q=bam_plp_auto(mplp, &tid, &pos, n_plp)))
        {
	  if (pos < beg || pos >= end) {continue;}
	  if(pmtr->out_file) {fprintf(ofp,"depth = %d\n",*n_plp);}
	  else {printf("depth = %d\n",*n_plp);}
	  j=1;
          for (j = 0; j < *n_plp; ++j){ 
          p = &q[j]; 
          if(pmtr->out_file) {
          fprintf(ofp,"%c",bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b),p->qpos)]);
          fprintf(ofp,"	%d",(bam1_qual(p->b)[p->qpos]));
          if((p->b->core.flag & BAM_FREVERSE)) {fprintf(ofp,"	-");} else  {fprintf(ofp,"	+");}
          if((p->b->core.flag & BAM_FREAD1)) {fprintf(ofp,"	rd1");} else  {fprintf(ofp,"	rd2");}
          fprintf(ofp,"	%d",(p->qpos+1));
          if((p->b->core.flag & BAM_FPROPER_PAIR)) {fprintf(ofp,"	pp");} else  {fprintf(ofp,"	npp");}
          uint8_t *tt = bam_aux_get(p->b, "RG");
          if(tt != NULL ) {tt++; fprintf(ofp,"	%s",tt);} else {fprintf(ofp,"	NA");}
          
          uint8_t *md = bam_aux_get(p->b, "MD");
          //md++;
	  if(md != NULL ) {md++;fprintf(ofp,"	%s", md);} else {fprintf(ofp,"	NA");}
          
          fprintf(ofp,"	%d",p->b->core.n_cigar);
          fprintf(ofp,"	%d",p->b->core.qual);
          if(p->b->core.isize < 0) {fprintf(ofp,"	%d",(p->b->core.isize)*(-1));} else {fprintf(ofp,"	%d",p->b->core.isize);} 
          fprintf(ofp,"	%s\n",bam1_qname(p->b));
          }
          else
          {
          printf("%c",bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b),p->qpos)]);
          printf("	%d",(bam1_qual(p->b)[p->qpos]));
          if((p->b->core.flag & BAM_FREVERSE)) {printf("	-");} else  {printf("	+");}
          if((p->b->core.flag & BAM_FREAD1)) {printf("	rd1");} else  {printf("	rd2");}
          printf("	%d",(p->qpos+1));
          if((p->b->core.flag & BAM_FPROPER_PAIR)) {printf("	pp");} else  {printf("	npp");}
          uint8_t *tt = bam_aux_get(p->b, "RG");
          //tt++;
          if(tt != NULL ) {tt++; printf("	%s",tt);} else {printf("	NA");}
          //printf("	%s",tt);
          uint8_t *md = bam_aux_get(p->b, "MD");
          //md++;
          if(md != NULL ) {md++;printf("	%s", md);} else {printf("	NA");}
          //printf("	%s", md);
          printf("	%d",p->b->core.n_cigar);
          printf("	%d",p->b->core.qual);
          if(p->b->core.isize < 0) {printf("	%d",(p->b->core.isize)*(-1));} else {printf("	%d",p->b->core.isize);} 
          printf("	%s\n",bam1_qname(p->b));
          }
          }   
       }
       if(!j) {if(pmtr->out_file) {fprintf(ofp,"depth = 0\n");} else {printf("depth = 0\n");}}
                }
                else
                {
                ifp=fopen(pmtr->pos_file,"r");
                while((ch=fgetc(ifp))!=EOF)
                {
                   j=0;
                   
                   tchr[j]=ch;j++;
                   while(((ch=fgetc(ifp))!='	'))
                   {tchr[j]=ch;j++;}
                   tchr[j]='\0';j=0;
                   
                   while(((ch=fgetc(ifp))!='	'))
                   {
                    if(ch != '\n')
                    {
                      tpos[j]=ch;j++;
                    }
                    else {break;}
                   }
                   
                   tpos[j]='\0';
                   if(ch != '\n') {while((ch=fgetc(ifp))!='\n');}
                   
                   sprintf(chrpos,"%s:%s-%s",tchr,tpos,tpos);
                
                   
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

        
        if(pmtr->out_file) { fprintf(ofp,"%s    ref_base = %c	",chrpos, toupper(ref[beg]));}
        else {printf("%s    ref_base = %c	",chrpos, toupper(ref[beg]));}
        n=0;
        if(pmtr->out_file) {
        for(i=10;i>0;i--) {fprintf(ofp,"%c",toupper(ref[beg-i])); if(toupper(ref[beg-i])=='C' || toupper(ref[beg-i])=='G') {n++;}}
        fprintf(ofp," %c ",toupper(ref[beg]));
        for(i=1;i<=10;i++) {fprintf(ofp,"%c",toupper(ref[beg+i])); if(toupper(ref[beg+i])=='C' || toupper(ref[beg+i])=='G') {n++;}}
        fprintf(ofp,"	GC%%[in 20bp flanks] = %.2f	", ((float)n/20.0)*100);
        }
        else
        {
        for(i=10;i>0;i--) {printf("%c",toupper(ref[beg-i])); if(toupper(ref[beg-i])=='C' || toupper(ref[beg-i])=='G') {n++;}}
        printf(" %c ",toupper(ref[beg]));
        for(i=1;i<=10;i++) {printf("%c",toupper(ref[beg+i])); if(toupper(ref[beg+i])=='C' || toupper(ref[beg+i])=='G') {n++;}}
        printf("	GC%%[in 20bp flanks] = %.2f	", ((float)n/20.0)*100);
        }
        
        
        
	mplp = bam_plp_init(read_bam, (void*)data); 
	n_plp = calloc(1, sizeof(int)); 

        j=0;
	while((q=bam_plp_auto(mplp, &tid, &pos, n_plp)))
        {
	  if (pos < beg || pos >= end) {continue;}
	  if(pmtr->out_file) {fprintf(ofp,"depth = %d\n",*n_plp);}
	  else {printf("depth = %d\n",*n_plp);}
	  j=1;
          for (j = 0; j < *n_plp; ++j){ 
          p = &q[j]; 
          if(pmtr->out_file) {
          fprintf(ofp,"%c",bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b),p->qpos)]);
          fprintf(ofp,"	%d",(bam1_qual(p->b)[p->qpos]));
          if((p->b->core.flag & BAM_FREVERSE)) {fprintf(ofp,"	-");} else  {fprintf(ofp,"	+");}
          if((p->b->core.flag & BAM_FREAD1)) {fprintf(ofp,"	rd1");} else  {fprintf(ofp,"	rd2");}
          fprintf(ofp,"	%d",(p->qpos+1));
          if((p->b->core.flag & BAM_FPROPER_PAIR)) {fprintf(ofp,"	pp");} else  {fprintf(ofp,"	npp");}
          uint8_t *tt = bam_aux_get(p->b, "RG");
          //tt++;
          if(tt != NULL ) {tt++;fprintf(ofp,"	%s",tt);} else {fprintf(ofp,"	NA");}
          
          uint8_t *md = bam_aux_get(p->b, "MD");
          //md++;
          if(md != NULL ) {md++;fprintf(ofp,"	%s", md);} else {fprintf(ofp,"	NA");}
          
          fprintf(ofp,"	%d",p->b->core.n_cigar);
          fprintf(ofp,"	%d",p->b->core.qual);
          if(p->b->core.isize < 0) {fprintf(ofp,"	%d",(p->b->core.isize)*(-1));} else {fprintf(ofp,"	%d",p->b->core.isize);} 
          fprintf(ofp,"	%s\n",bam1_qname(p->b));
          }
          else
          {
          printf("%c",bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b),p->qpos)]);
          printf("	%d",(bam1_qual(p->b)[p->qpos]));
          if((p->b->core.flag & BAM_FREVERSE)) {printf("	-");} else  {printf("	+");}
          if((p->b->core.flag & BAM_FREAD1)) {printf("	rd1");} else  {printf("	rd2");}
          printf("	%d",(p->qpos+1));
          if((p->b->core.flag & BAM_FPROPER_PAIR)) {printf("	pp");} else  {printf("	npp");}
          uint8_t *tt = bam_aux_get(p->b, "RG");
          //tt++;
          if(tt != NULL ) {tt++;printf("	%s",tt);} else {printf("	NA");}
          
          uint8_t *md = bam_aux_get(p->b, "MD");
          //md++;
          if(md != NULL ) {md++; printf("	%s", md);} else {printf("	NA");}
         // printf("	%s", md);
          printf("	%d",p->b->core.n_cigar);
          printf("	%d",p->b->core.qual);
          if(p->b->core.isize < 0) {printf("	%d",(p->b->core.isize)*(-1));} else {printf("	%d",p->b->core.isize);} 
          printf("	%s\n",bam1_qname(p->b));
          }
          } 
        }
       if(!j) {if(pmtr->out_file) {fprintf(ofp,"depth = 0\n");} else {printf("depth = 0\n");}}	
}
fclose(ifp);
}		
if(pmtr->out_file) {fclose(ofp);}			
	
		
	
	
	bam_plp_destroy(mplp);
        bam_index_destroy(idx);
	bam_header_destroy(h);
	bam_close(data->fp);
	bam_iter_destroy(data->iter);
	free(ref);
	free(data); 
	return 0;
}

