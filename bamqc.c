
 #include <stdio.h>
 #include <stdlib.h>
 #include <time.h>
 #include "bam.h"
 #include "faidx.h"
 #include "kstring.h"

 #include "ksort.h"
KSORT_INIT_GENERIC(uint64_t)

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 8192)

 typedef struct {
    int n, m;
    uint64_t *a;
    int *idx;
    int filter;
} bed_reglist_t;
 #include "khash.h"
KHASH_MAP_INIT_STR(reg, bed_reglist_t)
#define LIDX_SHIFT 13
typedef kh_reg_t reghash_t;

 #define _cop(c) ((c)&BAM_CIGAR_MASK)
 #define _cln(c) ((c)>>BAM_CIGAR_SHIFT)
 
 #define AT 0
 #define AC 1
 #define AG 2
 #define TA 3
 #define TC 4
 #define TG 5
 #define CA 6
 #define CT 7
 #define CG 8
 #define GA 9
 #define GT 10
 #define GC 11

 #define AA 12
 #define TT 13
 #define GG 14
 #define CC 15


 #define BAM_FSUPPLEMENTARY 2048
 
 #define MAXREADLENGTH 250
 #define MAX_ISIZE 600
 #define MAX_BASEQ 40
 #define MAX_MAPQ 60
 
 char Adp_prim_trans[8][14] = {"CTGTCTCTTATA" , "CAAGCAGAAGAC" , "GTCTTCTGCTTG" , "AGATCGGAAGAG" , "CTCTTCCGATCT" , "AATGATACGGCG" };
 
typedef struct {
    char *outdir, *ref_file, *bed_file, *bam_file, *xrg_file;
    char *sample_name, *comnd, *library;
    char sId[500][1024], x_rg[500][1024];
    long int ref_size, trgtsize;
    int thrd, no_of_contigs, no_of_rgs, xrg_cnt;
} paramtr;


typedef struct {
    int dif_count;
    int q_sum;
    int sc_pct;
    int sq_len;
} read_info_t;

typedef struct {
long int rd_posi[MAXREADLENGTH];
} RP;

typedef struct {
int A,T,G,C;
} ATGC;

typedef struct {
long int basemmc;
long int polyIns, polyDel;
long int del_cnt, del_cnt_bs;
long int ins_cnt, ins_cnt_bs;
long int hrdclp, hrdclp_bs;
long int sftclp, sftclp_bs;
} MIDCLP;

typedef struct {
long int NTS, IPCRprm, IADptr;
} ADPCR;

typedef struct {
   long int A_con[MAXREADLENGTH], T_con[MAXREADLENGTH], G_con[MAXREADLENGTH], C_con[MAXREADLENGTH], N_con[MAXREADLENGTH], Percycle_mn_baseQ[MAXREADLENGTH];
   long int r1_cnt, r2_cnt, rev_cnt, map_pr, propr_pr, qc_fail;
   long int supply_cnt, scndry_cnt, unmap, all_rd_cnt, tot_rd_cnt, map_difchr, DupCnt;
   long int GC_pct[101], BaseQ_dist[MAX_BASEQ+1],Isize_dist[MAX_ISIZE],MapQ_dist[MAX_MAPQ], MapQ_bin_dist[10], Rdlen[MAXREADLENGTH];
   long int AAcnt, TTcnt, GGcnt, CCcnt, NNcnt, FF, FR, RR, OO, Isize_ge1K;
} log_genrl;

typedef struct {
int qv1,qv2,qv3,qv8;
} QV_range;

typedef struct {
    int tmmc[6];
    int tDel[5],tIns[5];
    int tDel_flnks[5][16];
    int tIns_flnks[5][16];
    int tClp3[5],tClp5[5];
    int MMC[12];
    int MMC_flnks[12][16];
    ATGC Del1,Ins1;
    QV_range mmqv[12];
    RP rmmc[12], rmmQV;
    MIDCLP midc;
    ADPCR adptr;
    long int no_indclp;
    long long int cnt;
    log_genrl glog;
} log_v;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
int *bed_index_core_u(int n, uint64_t *a, int *n_idx)
{
        int i, j, m, *idx;
        m = *n_idx = 0; idx = 0;
        for (i = 0; i < n; ++i) {
                int beg, end;
                beg = a[i]>>32 >> LIDX_SHIFT; end = ((uint32_t)a[i]) >> LIDX_SHIFT;
                if (m < end + 1) {
                        int oldm = m;
                        m = end + 1;
                        kroundup32(m);
                        idx = realloc(idx, m * sizeof(int));
                        for (j = oldm; j < m; ++j) idx[j] = -1;
                }
                if (beg == end) {
                        if (idx[beg] < 0) idx[beg] = i;
                } else {
                        for (j = beg; j <= end; ++j)
                                if (idx[j] < 0) idx[j] = i;
                }
                *n_idx = end + 1;
        }
        return idx;
}

void bed_index_u(void *_h)
{
        reghash_t *h = (reghash_t*)_h;
        khint_t k;
        for (k = 0; k < kh_end(h); ++k) {
                if (kh_exist(h, k)) {
                        bed_reglist_t *p = &kh_val(h, k);
                        if (p->idx) free(p->idx);
                        ks_introsort(uint64_t, p->n, p->a);
                        p->idx = bed_index_core_u(p->n, p->a, &p->m);
                }
        }
}

int bed_overlap_core_u(const bed_reglist_t *p, int beg, int end, int *read_ovrlp_st, int *read_ovrlp_nd)
{
	int i, min_off;
	//printf("\nrdst=%ld rdnd=%ld\n", *read_ovrlp_st, *read_ovrlp_nd);
	if (p->n == 0) return 0;
	min_off = (beg>>LIDX_SHIFT >= p->n)? p->idx[p->n-1] : p->idx[beg>>LIDX_SHIFT];
	if (min_off < 0) { // TODO: this block can be improved, but speed should not matter too much here
		int n = beg>>LIDX_SHIFT;
		if (n > p->n) n = p->n;
		for (i = n - 1; i >= 0; --i)
			if (p->idx[i] >= 0) break;
		min_off = i >= 0? p->idx[i] : 0;
	}
	for (i = min_off; i < p->n; ++i) {
		if ((int)(p->a[i]>>32) >= end) break; // out of range; no need to proceed
		if ((int32_t)p->a[i] > beg && (int32_t)(p->a[i]>>32) < end){
                if (beg < (int32_t)(p->a[i]>>32)) 
                { 
                    *read_ovrlp_st = (int32_t)(p->a[i]>>32) - beg; /* bed_start - read_start_pos*/
                }
                if (end > (int32_t)p->a[i])
                {
                    //*read_ovrlp_st = 0; 
                    //printf("I am in end condition\n");
                    *read_ovrlp_nd = (((int32_t)(p->a[i])) - beg); /*bed_end - read_end_pos*/
                    //*read_ovrlp_nd = end - (int32_t)p->a[i];
                }
                
               //printf("bdst = %d	bdnd = %d	mapst=%d	mapnd=%d\n", (int32_t)(p->a[i]>>32), (int32_t)p->a[i], beg, end);
                        
			return 1;} // find the overlap; return
	}
	return 0;
}

int bed_overlap_u(const void *_h, const char *chr, int beg, int end, int *read_ovrlp_st, int *read_ovrlp_nd)
{
	const reghash_t *h = (const reghash_t*)_h;
	khint_t k;
	if (!h) return 0;
	k = kh_get(reg, h, chr);
	if (k == kh_end(h)) return 0;
	return bed_overlap_core_u(&kh_val(h, k), beg, end, read_ovrlp_st, read_ovrlp_nd);
}

void *bed_read_new(const char *fn, long int *trgtsize)
{
        reghash_t *h = kh_init(reg);
        gzFile fp;
        kstream_t *ks;
        uint64_t sum=0;
        int dret;
        kstring_t *str;
        //*trgtsize = 5;
        fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
        if (fp == 0) { printf("Unable to open bed file\n"); return 0;}
        str = calloc(1, sizeof(kstring_t));
        ks = ks_init(fp);
        while (ks_getuntil(ks, 0, str, &dret) >= 0) { 
                int beg = -1, end = -1;
                bed_reglist_t *p;
                khint_t k = kh_get(reg, h, str->s);
                if (k == kh_end(h)) { 
                        int ret;
                        char *s = strdup(str->s);
                        k = kh_put(reg, h, s, &ret);
                        memset(&kh_val(h, k), 0, sizeof(bed_reglist_t));
                }
                p = &kh_val(h, k);
                //printf("p->n=%d	p->m=%d\n",p->n,p->m);
                if (dret != '\n') { 
                        if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
                                beg = atoi(str->s); 
                                if (dret != '\n') {
                                        if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
                                                end = atoi(str->s); 
                                                if (end < beg) end = -1;
                                        }
                                }
                        }
                }
                if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n'); 
                if (end < 0 && beg > 0) end = beg, beg = beg - 1; 
                if (beg >= 0 && end > beg) {
                        sum+=(end-beg);
                        if (p->n == p->m) {
                                p->m = p->m? p->m<<1 : 4;
                                p->a = realloc(p->a, p->m * 8);
                        }
                        p->a[p->n++] = (uint64_t)beg<<32 | end;
                }
        }
        *trgtsize = sum;
        ks_destroy(ks);
        gzclose(fp);
        free(str->s); free(str);
        bed_index_u(h);
        return h;
}
void bed_destroy_u(void *_h)
{
        reghash_t *h = (reghash_t*)_h;
        khint_t k;
        for (k = 0; k < kh_end(h); ++k) {
                if (kh_exist(h, k)) {
                        free(kh_val(h, k).a);
                        free(kh_val(h, k).idx);
                        free((char*)kh_key(h, k));
                }
        }
        kh_destroy(reg, h);
}

static void bed_print1(void *reg_hash) {
    reghash_t *h = (reghash_t *)reg_hash;
    bed_reglist_t *p;
    khint_t k;
    int i;
    const char *reg;
    uint32_t beg, end;

    if (!h) {
        printf("Hash table is empty!\n");
        return;
    }
    for (k = kh_begin(h); k < kh_end(h); k++) {
        if (kh_exist(h,k)) {
            reg = kh_key(h,k);
            //printf("Region: '%s'\n", reg);
            p = &kh_val(h,k);
            //printf("p->n=%d	p->m=%d\n",p->n,p->m);
            if ((p = &kh_val(h,k)) != NULL && p->n > 0) {
                //printf("Filter: %d\n", p->filter);
                for (i=0; i<(p->n); i++) {
                    beg = (uint32_t)(p->a[i]>>32);
                    end = (uint32_t)(p->a[i]);

                    printf("%s	%d	%d\n",reg,beg,end);
                }
            } else {
                printf("Region '%s' has no intervals!\n", reg);
            }
        }
    }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//




void iniT(log_v *t)
{
  int i,j;
  for(i=0;i<5;i++) 
  {
      t->tmmc[i]=0; t->tDel[i]=0; t->tIns[i]=0;
      t->tClp3[i]=0; t->tClp5[i]=0;
  }
  t->Del1.A=0; t->Del1.T=0;t->Del1.G=0; t->Del1.C=0;
  t->Ins1.A=0; t->Ins1.T=0;t->Ins1.G=0; t->Ins1.C=0;
  for(i=0;i<12;i++) {t->mmqv[i].qv1=0; t->mmqv[i].qv2=0; t->mmqv[i].qv3=0; t->mmqv[i].qv8=0;}
  for(i=0;i<12;i++) {for(j=0;j<250;j++){t->rmmc[i].rd_posi[j]=0;}}
}

void Addup_u_u(log_v *OA, log_v *r12, log_v *rfr, log_v *t, int flag, int rLen, int revorfor/*0 for forward and 1 for reverse*/)
{
  int i,j;
  //OA->cnt+=t->cnt; r12->cnt+=t->cnt; rfr->cnt+=t->cnt;
  OA->no_indclp+=t->no_indclp; r12->no_indclp+=t->no_indclp; rfr->no_indclp+=t->no_indclp;
  OA->glog.tot_rd_cnt+=t->glog.tot_rd_cnt; r12->glog.tot_rd_cnt+=t->glog.tot_rd_cnt; rfr->glog.tot_rd_cnt+=t->glog.tot_rd_cnt;
  
  if(flag) {t->cnt=0;t->no_indclp=0; t->glog.tot_rd_cnt=0;}
  
  OA->midc.basemmc+=t->midc.basemmc; r12->midc.basemmc+=t->midc.basemmc; // rtr->midc.basemmc+=t->midc.basemmc;
  OA->midc.polyIns+=t->midc.polyIns; r12->midc.polyIns+=t->midc.polyIns; // rfr->midc.polyIns+=t->midc.polyIns; 
  OA->midc.polyDel+=t->midc.polyDel; r12->midc.polyDel+=t->midc.polyDel; // rfr->midc.polyDel+=t->midc.polyDel;
  OA->midc.del_cnt+=t->midc.del_cnt; r12->midc.del_cnt+=t->midc.del_cnt; // rfr->midc.del_cnt+=t->midc.del_cnt; 
  OA->midc.del_cnt_bs+=t->midc.del_cnt_bs; r12->midc.del_cnt_bs+=t->midc.del_cnt_bs; // rfr->midc.del_cnt_bs+=t->midc.del_cnt_bs;
  OA->midc.ins_cnt+=t->midc.ins_cnt; r12->midc.ins_cnt+=t->midc.ins_cnt; // rfr->midc.ins_cnt+=t->midc.ins_cnt; 
  OA->midc.ins_cnt_bs+=t->midc.ins_cnt_bs; r12->midc.ins_cnt_bs+=t->midc.ins_cnt_bs; // rfr->midc.ins_cnt_bs+=t->midc.ins_cnt_bs;
  OA->midc.hrdclp+=t->midc.hrdclp; r12->midc.hrdclp+=t->midc.hrdclp; // rfr->midc.hrdclp+=t->midc.hrdclp; 
  OA->midc.hrdclp_bs+=t->midc.hrdclp_bs; r12->midc.hrdclp_bs+=t->midc.hrdclp_bs; // rfr->midc.hrdclp_bs+=t->midc.hrdclp_bs;
  OA->midc.sftclp+=t->midc.sftclp; r12->midc.sftclp+=t->midc.sftclp; // rfr->midc.sftclp+=t->midc.sftclp; 
  OA->midc.sftclp_bs+=t->midc.sftclp_bs; r12->midc.sftclp_bs+=t->midc.sftclp_bs; // rfr->midc.sftclp_bs+=t->midc.sftclp_bs;
  
  OA->adptr.IADptr+=t->adptr.IADptr; r12->adptr.IADptr+=t->adptr.IADptr; // rfr->adptr.IADptr+=t->adptr.IADptr; 
  OA->adptr.IPCRprm+=t->adptr.IPCRprm; r12->adptr.IPCRprm+=t->adptr.IPCRprm; // rfr->adptr.IPCRprm+=t->adptr.IPCRprm; 
  OA->adptr.NTS+=t->adptr.NTS; r12->adptr.NTS+=t->adptr.NTS; // rfr->adptr.NTS+=t->adptr.NTS;
  
  OA->glog.r1_cnt+=t->glog.r1_cnt; r12->glog.r1_cnt+=t->glog.r1_cnt; // rfr->glog.r1_cnt+=t->glog.r1_cnt; 
  OA->glog.r2_cnt+=t->glog.r2_cnt; r12->glog.r2_cnt+=t->glog.r2_cnt; // rfr->glog.r2_cnt+=t->glog.r2_cnt; 
  OA->glog.rev_cnt+=t->glog.rev_cnt; r12->glog.rev_cnt+=t->glog.rev_cnt; // rfr->glog.rev_cnt+=t->glog.rev_cnt;
  OA->glog.map_pr+=t->glog.map_pr; r12->glog.map_pr+=t->glog.map_pr; // rfr->glog.map_pr+=t->glog.map_pr; 
  OA->glog.propr_pr+=t->glog.propr_pr; r12->glog.propr_pr+=t->glog.propr_pr; // rfr->glog.propr_pr+=t->glog.propr_pr; 
  OA->glog.qc_fail+=t->glog.qc_fail; r12->glog.qc_fail+=t->glog.qc_fail; // rfr->glog.qc_fail+=t->glog.qc_fail;
  OA->glog.supply_cnt+=t->glog.supply_cnt; r12->glog.supply_cnt+=t->glog.supply_cnt; // rfr->glog.supply_cnt+=t->glog.supply_cnt; 
  OA->glog.scndry_cnt+=t->glog.scndry_cnt; r12->glog.scndry_cnt+=t->glog.scndry_cnt; // rfr->glog.scndry_cnt+=t->glog.scndry_cnt; 
  
  OA->glog.AAcnt+=t->glog.AAcnt; r12->glog.AAcnt+=t->glog.AAcnt; // rfr->glog.AAcnt+=t->glog.AAcnt; 
  OA->glog.TTcnt+=t->glog.TTcnt; r12->glog.TTcnt+=t->glog.TTcnt; // rfr->glog.TTcnt+=t->glog.TTcnt; 
  OA->glog.GGcnt+=t->glog.GGcnt; r12->glog.GGcnt+=t->glog.GGcnt; // rfr->glog.GGcnt+=t->glog.GGcnt;
  OA->glog.CCcnt+=t->glog.CCcnt; r12->glog.CCcnt+=t->glog.CCcnt; // rfr->glog.CCcnt+=t->glog.CCcnt; 
  OA->glog.NNcnt+=t->glog.NNcnt; r12->glog.NNcnt+=t->glog.NNcnt; // rfr->glog.NNcnt+=t->glog.NNcnt;
  
  OA->glog.map_difchr+=t->glog.map_difchr;
  OA->glog.DupCnt+=t->glog.DupCnt;
  OA->glog.Isize_ge1K+=t->glog.Isize_ge1K;
  
  OA->glog.FF+=t->glog.FF; 
  OA->glog.FR+=t->glog.FR;
  OA->glog.RR+=t->glog.RR; 
  OA->glog.OO+=t->glog.OO;
  
  if(flag)
  {
    t->midc.basemmc=0;
    t->midc.polyIns=0; t->midc.polyDel=0;
    t->midc.del_cnt=0; t->midc.del_cnt_bs=0;
    t->midc.ins_cnt=0; t->midc.ins_cnt_bs=0;
    t->midc.hrdclp=0;  t->midc.hrdclp_bs=0;
    t->midc.sftclp=0;  t->midc.sftclp_bs=0;
  
    t->adptr.IADptr=0; t->adptr.IPCRprm=0; t->adptr.NTS=0;
  
    t->glog.r1_cnt=0; t->glog.r2_cnt=0; t->glog.rev_cnt=0;
    t->glog.map_pr=0; t->glog.propr_pr=0; t->glog.qc_fail=0;
    t->glog.supply_cnt=0; t->glog.scndry_cnt=0; t->glog.map_difchr=0;
    t->glog.DupCnt=0; t->glog.Isize_ge1K=0;
  
    t->glog.AAcnt=0; t->glog.TTcnt=0; t->glog.GGcnt=0;
    t->glog.CCcnt=0; t->glog.NNcnt=0;
  
    t->glog.FF=0; t->glog.FR=0;
    t->glog.RR=0; t->glog.OO=0;
  }
  
  for(i=0;i<rLen;i++)
  {
    OA->glog.A_con[i]+=t->glog.A_con[i];
    OA->glog.T_con[i]+=t->glog.T_con[i];
    OA->glog.G_con[i]+=t->glog.G_con[i];
    OA->glog.C_con[i]+=t->glog.C_con[i];
    OA->glog.N_con[i]+=t->glog.N_con[i];
    OA->glog.Rdlen[i]+=t->glog.Rdlen[i];
    OA->glog.Percycle_mn_baseQ[i]+=t->glog.Percycle_mn_baseQ[i];
    
    if(flag)
    {
       t->glog.A_con[i]=0;
       t->glog.T_con[i]=0;
       t->glog.G_con[i]=0;
       t->glog.C_con[i]=0;
       t->glog.N_con[i]=0;
       t->glog.Rdlen[i]=0;
       t->glog.Percycle_mn_baseQ[i]=0;
    }
  }

  
  for(i=0;i<5;i++) 
  {
      OA->tmmc[i]+=t->tmmc[i]; rfr->tmmc[i]+=t->tmmc[i];
      OA->tDel[i]+=t->tDel[i]; rfr->tDel[i]+=t->tDel[i];
      OA->tIns[i]+=t->tIns[i]; rfr->tIns[i]+=t->tIns[i];
      OA->tClp3[i]+=t->tClp3[i]; rfr->tClp3[i]+=t->tClp3[i];
      OA->tClp5[i]+=t->tClp5[i]; rfr->tClp5[i]+=t->tClp5[i];
      
      for(j=0;j<16;j++)
      {
         OA->tDel_flnks[i][j]+=t->tDel_flnks[i][j];
         OA->tIns_flnks[i][j]+=t->tIns_flnks[i][j];
         rfr->tDel_flnks[i][j]+=t->tDel_flnks[i][j];
         rfr->tIns_flnks[i][j]+=t->tIns_flnks[i][j];
         if(flag)
         {
           t->tDel_flnks[i][j]=0; t->tIns_flnks[i][j]=0;
         }
      }
      
      if(flag)
      {
         t->tmmc[i]=0;
         t->tDel[i]=0;
         t->tIns[i]=0;
         t->tClp3[i]=0;
         t->tClp5[i]=0;
      }
  }
  OA->Del1.A+=t->Del1.A; OA->Del1.T+=t->Del1.T; OA->Del1.G+=t->Del1.G; OA->Del1.C+=t->Del1.C;
  OA->Ins1.A+=t->Ins1.A; OA->Ins1.T+=t->Ins1.T; OA->Ins1.G+=t->Ins1.G; OA->Ins1.C+=t->Ins1.C;

  rfr->Del1.A+=t->Del1.A; rfr->Del1.T+=t->Del1.T; rfr->Del1.G+=t->Del1.G; rfr->Del1.C+=t->Del1.C;
  rfr->Ins1.A+=t->Ins1.A; rfr->Ins1.T+=t->Ins1.T; rfr->Ins1.G+=t->Ins1.G; rfr->Ins1.C+=t->Ins1.C;
  if(flag)
  {
      t->Del1.A=0; t->Del1.T=0; t->Del1.G=0; t->Del1.C=0;
      t->Ins1.A=0; t->Ins1.T=0; t->Ins1.G=0; t->Ins1.C=0;
  }
  for(i=0;i<12;i++) 
  {
      OA->mmqv[i].qv1+=t->mmqv[i].qv1; OA->mmqv[i].qv2+=t->mmqv[i].qv2;
      OA->mmqv[i].qv3+=t->mmqv[i].qv3; OA->mmqv[i].qv8+=t->mmqv[i].qv8;

      rfr->mmqv[i].qv1+=t->mmqv[i].qv1; rfr->mmqv[i].qv2+=t->mmqv[i].qv2;
      rfr->mmqv[i].qv3+=t->mmqv[i].qv3; rfr->mmqv[i].qv8+=t->mmqv[i].qv8;

      if(flag)
      {
         t->mmqv[i].qv1=0; t->mmqv[i].qv2=0;
         t->mmqv[i].qv3=0; t->mmqv[i].qv8=0;
      }

  }

  for(i=0;i<12;i++) 
  {
      OA->MMC[i]+=t->MMC[i];
      rfr->MMC[i]+=t->MMC[i];
      if(flag) { t->MMC[i]=0; }
  } 
  
  for(i=0;i<12;i++) 
  {
     for(j=0;j<16;j++)
     {
      OA->MMC_flnks[i][j]+=t->MMC_flnks[i][j];
      rfr->MMC_flnks[i][j]+=t->MMC_flnks[i][j];
      if(flag) { t->MMC_flnks[i][j]=0; }
     }
  }
  
  
   
 if(revorfor == 0) // if forward
 {
  for(j=0;j<rLen;j++)
  {
      OA->rmmQV.rd_posi[j]+=t->rmmQV.rd_posi[j];
      rfr->rmmQV.rd_posi[j]+=t->rmmQV.rd_posi[j];
      if(flag) {t->rmmQV.rd_posi[j]=0;}
  }
  
  for(i=0;i<12;i++)
  {
      for(j=0;j<rLen;j++)
      {
          OA->rmmc[i].rd_posi[j]+=t->rmmc[i].rd_posi[j]; // Total substitution rate
          rfr->rmmc[i].rd_posi[j]+=t->rmmc[i].rd_posi[j]; 
          if(flag) {t->rmmc[i].rd_posi[j]=0;}
      }
  }
 }
 else // if reverse
 {
  for(j=0;j<rLen;j++)
  {
      OA->rmmQV.rd_posi[j]+=t->rmmQV.rd_posi[(rLen-1)-j];
      rfr->rmmQV.rd_posi[j]+=t->rmmQV.rd_posi[(rLen-1)-j];
      if(flag) {t->rmmQV.rd_posi[(rLen-1)-j]=0;}
  }
  
  for(i=0;i<12;i++)
  {
      for(j=0;j<rLen;j++)
      {
          OA->rmmc[i].rd_posi[j]+=t->rmmc[i].rd_posi[(rLen-1)-j]; // Total substitution rate
          rfr->rmmc[i].rd_posi[j]+=t->rmmc[i].rd_posi[(rLen-1)-j]; 
          if(flag) {t->rmmc[i].rd_posi[(rLen-1)-j]=0;}
      }
  }
 } 
  
}

void Addup_u(log_v *d, log_v *t, int flag, int rLen, char *prefix)
{
  int i,j;
  d->cnt+=t->cnt;
  d->no_indclp+=t->no_indclp;
  d->glog.tot_rd_cnt+=t->glog.tot_rd_cnt;
  if(flag) {t->cnt=0;t->no_indclp=0; t->glog.tot_rd_cnt=0;}
  
  d->midc.basemmc+=t->midc.basemmc;
  d->midc.polyIns+=t->midc.polyIns; d->midc.polyDel+=t->midc.polyDel;
  d->midc.del_cnt+=t->midc.del_cnt; d->midc.del_cnt_bs+=t->midc.del_cnt_bs;
  d->midc.ins_cnt+=t->midc.ins_cnt; d->midc.ins_cnt_bs+=t->midc.ins_cnt_bs;
  d->midc.hrdclp+=t->midc.hrdclp; d->midc.hrdclp_bs+=t->midc.hrdclp_bs;
  d->midc.sftclp+=t->midc.sftclp; d->midc.sftclp_bs+=t->midc.sftclp_bs;
  
  d->adptr.IADptr+=t->adptr.IADptr; d->adptr.IPCRprm+=t->adptr.IPCRprm; d->adptr.NTS+=t->adptr.NTS;
  
  d->glog.r1_cnt+=t->glog.r1_cnt; d->glog.r2_cnt+=t->glog.r2_cnt; d->glog.rev_cnt+=t->glog.rev_cnt;
  d->glog.map_pr+=t->glog.map_pr; d->glog.propr_pr+=t->glog.propr_pr; d->glog.qc_fail+=t->glog.qc_fail;
  d->glog.supply_cnt+=t->glog.supply_cnt; d->glog.scndry_cnt+=t->glog.scndry_cnt; d->glog.map_difchr+=t->glog.map_difchr;
  d->glog.DupCnt+=t->glog.DupCnt; d->glog.Isize_ge1K+=t->glog.Isize_ge1K;
  
  d->glog.AAcnt+=t->glog.AAcnt; d->glog.TTcnt+=t->glog.TTcnt; d->glog.GGcnt+=t->glog.GGcnt;
  d->glog.CCcnt+=t->glog.CCcnt; d->glog.NNcnt+=t->glog.NNcnt;
  
  d->glog.FF+=t->glog.FF; d->glog.FR+=t->glog.FR;
  d->glog.RR+=t->glog.RR; d->glog.OO+=t->glog.OO;
  
  if(flag)
  {
    t->midc.basemmc=0;
    t->midc.polyIns=0; t->midc.polyDel=0;
    t->midc.del_cnt=0; t->midc.del_cnt_bs=0;
    t->midc.ins_cnt=0; t->midc.ins_cnt_bs=0;
    t->midc.hrdclp=0;  t->midc.hrdclp_bs=0;
    t->midc.sftclp=0;  t->midc.sftclp_bs=0;
  
    t->adptr.IADptr=0; t->adptr.IPCRprm=0; t->adptr.NTS=0;
  
    t->glog.r1_cnt=0; t->glog.r2_cnt=0; t->glog.rev_cnt=0;
    t->glog.map_pr=0; t->glog.propr_pr=0; t->glog.qc_fail=0;
    t->glog.supply_cnt=0; t->glog.scndry_cnt=0; t->glog.map_difchr=0;
    t->glog.DupCnt=0; t->glog.Isize_ge1K=0;
  
    t->glog.AAcnt=0; t->glog.TTcnt=0; t->glog.GGcnt=0;
    t->glog.CCcnt=0; t->glog.NNcnt=0;
  
    t->glog.FF=0; t->glog.FR=0;
    t->glog.RR=0; t->glog.OO=0;
  }
  
  for(i=0;i<rLen;i++)
  {
    d->glog.A_con[i]+=t->glog.A_con[i];
    d->glog.T_con[i]+=t->glog.T_con[i];
    d->glog.G_con[i]+=t->glog.G_con[i];
    d->glog.C_con[i]+=t->glog.C_con[i];
    d->glog.N_con[i]+=t->glog.N_con[i];
    d->glog.Rdlen[i]+=t->glog.Rdlen[i];
    d->glog.Percycle_mn_baseQ[i]+=t->glog.Percycle_mn_baseQ[i];
    
    if(flag)
    {
       t->glog.A_con[i]=0;
       t->glog.T_con[i]=0;
       t->glog.G_con[i]=0;
       t->glog.C_con[i]=0;
       t->glog.N_con[i]=0;
       t->glog.Rdlen[i]=0;
       t->glog.Percycle_mn_baseQ[i]=0;
    }
  }

  
  for(i=0;i<5;i++) 
  {
      d->tmmc[i]+=t->tmmc[i];
      d->tDel[i]+=t->tDel[i];
      d->tIns[i]+=t->tIns[i];
      d->tClp3[i]+=t->tClp3[i];
      d->tClp5[i]+=t->tClp5[i];
      
      if(flag)
      {
         t->tmmc[i]=0;
         t->tDel[i]=0;
         t->tIns[i]=0;
         t->tClp3[i]=0;
         t->tClp5[i]=0;
      }
  }
  d->Del1.A+=t->Del1.A; d->Del1.T+=t->Del1.T; d->Del1.G+=t->Del1.G; d->Del1.C+=t->Del1.C;
  d->Ins1.A+=t->Ins1.A; d->Ins1.T+=t->Ins1.T; d->Ins1.G+=t->Ins1.G; d->Ins1.C+=t->Ins1.C;
  if(flag)
  {
      t->Del1.A=0; t->Del1.T=0; t->Del1.G=0; t->Del1.C=0;
      t->Ins1.A=0; t->Ins1.T=0; t->Ins1.G=0; t->Ins1.C=0;
  }
  for(i=0;i<12;i++) 
  {
      d->mmqv[i].qv1+=t->mmqv[i].qv1; d->mmqv[i].qv2+=t->mmqv[i].qv2;
      d->mmqv[i].qv3+=t->mmqv[i].qv3; d->mmqv[i].qv8+=t->mmqv[i].qv8;

      if(flag)
      {
         t->mmqv[i].qv1=0; t->mmqv[i].qv2=0;
         t->mmqv[i].qv3=0; t->mmqv[i].qv8=0;
      }

  }

  for(i=0;i<12;i++) 
  {
      d->MMC[i]+=t->MMC[i];

      if(flag) { t->MMC[i]=0; }
  }  
  

  for(j=0;j<rLen;j++)
  {
      d->rmmQV.rd_posi[j]+=t->rmmQV.rd_posi[j];
      if(flag) {t->rmmQV.rd_posi[j]=0;}
  }
  
  for(i=0;i<12;i++)
  {
      for(j=0;j<rLen;j++)
      {
          d->rmmc[i].rd_posi[j]+=t->rmmc[i].rd_posi[j];
          if(flag) {t->rmmc[i].rd_posi[j]=0;}
      }
  }
}
void binQual(uint8_t qv, int index, log_v *T)
{
  if(qv <= 10) {T->mmqv[index].qv1+=1;}
  else if(qv > 10 && qv <= 20) {T->mmqv[index].qv2+=1;}
  else if(qv > 20 && qv <= 30) {T->mmqv[index].qv3+=1;}
  else if(qv > 30) {T->mmqv[index].qv8+=1;}
}

void count_flanks(char pre, char post, int chang, int flag, log_v *T)
{
  if(flag == 0)
  {
       if(pre=='A' && post=='A') {T->MMC_flnks[chang][AA]+=1;}
  else if(pre=='A' && post=='T') {T->MMC_flnks[chang][AT]+=1;}
  else if(pre=='A' && post=='G') {T->MMC_flnks[chang][AG]+=1;}
  else if(pre=='A' && post=='C') {T->MMC_flnks[chang][AC]+=1;}
  else if(pre=='T' && post=='A') {T->MMC_flnks[chang][TA]+=1;}
  else if(pre=='T' && post=='T') {T->MMC_flnks[chang][TT]+=1;}
  else if(pre=='T' && post=='G') {T->MMC_flnks[chang][TG]+=1;}
  else if(pre=='T' && post=='C') {T->MMC_flnks[chang][TC]+=1;}
  else if(pre=='G' && post=='A') {T->MMC_flnks[chang][GA]+=1;}
  else if(pre=='G' && post=='T') {T->MMC_flnks[chang][GT]+=1;}
  else if(pre=='G' && post=='G') {T->MMC_flnks[chang][GG]+=1;}
  else if(pre=='G' && post=='C') {T->MMC_flnks[chang][GC]+=1;}
  else if(pre=='C' && post=='A') {T->MMC_flnks[chang][CA]+=1;}
  else if(pre=='C' && post=='T') {T->MMC_flnks[chang][CT]+=1;}
  else if(pre=='C' && post=='G') {T->MMC_flnks[chang][CG]+=1;}
  else if(pre=='C' && post=='C') {T->MMC_flnks[chang][CC]+=1;}
  }
  else if(flag == 1)
  { 
       if(pre=='A' && post=='A') {T->tDel_flnks[chang][AA]+=1;}
  else if(pre=='A' && post=='T') {T->tDel_flnks[chang][AT]+=1;}
  else if(pre=='A' && post=='G') {T->tDel_flnks[chang][AG]+=1;}
  else if(pre=='A' && post=='C') {T->tDel_flnks[chang][AC]+=1;}
  else if(pre=='T' && post=='A') {T->tDel_flnks[chang][TA]+=1;}
  else if(pre=='T' && post=='T') {T->tDel_flnks[chang][TT]+=1;}
  else if(pre=='T' && post=='G') {T->tDel_flnks[chang][TG]+=1;}
  else if(pre=='T' && post=='C') {T->tDel_flnks[chang][TC]+=1;}
  else if(pre=='G' && post=='A') {T->tDel_flnks[chang][GA]+=1;}
  else if(pre=='G' && post=='T') {T->tDel_flnks[chang][GT]+=1;}
  else if(pre=='G' && post=='G') {T->tDel_flnks[chang][GG]+=1;}
  else if(pre=='G' && post=='C') {T->tDel_flnks[chang][GC]+=1;}
  else if(pre=='C' && post=='A') {T->tDel_flnks[chang][CA]+=1;}
  else if(pre=='C' && post=='T') {T->tDel_flnks[chang][CT]+=1;}
  else if(pre=='C' && post=='G') {T->tDel_flnks[chang][CG]+=1;}
  else if(pre=='C' && post=='C') {T->tDel_flnks[chang][CC]+=1;}
  }
  else if(flag == 2)
  {
       if(pre=='A' && post=='A') {T->tIns_flnks[chang][AA]+=1;}
  else if(pre=='A' && post=='T') {T->tIns_flnks[chang][AT]+=1;}
  else if(pre=='A' && post=='G') {T->tIns_flnks[chang][AG]+=1;}
  else if(pre=='A' && post=='C') {T->tIns_flnks[chang][AC]+=1;}
  else if(pre=='T' && post=='A') {T->tIns_flnks[chang][TA]+=1;}
  else if(pre=='T' && post=='T') {T->tIns_flnks[chang][TT]+=1;}
  else if(pre=='T' && post=='G') {T->tIns_flnks[chang][TG]+=1;}
  else if(pre=='T' && post=='C') {T->tIns_flnks[chang][TC]+=1;}
  else if(pre=='G' && post=='A') {T->tIns_flnks[chang][GA]+=1;}
  else if(pre=='G' && post=='T') {T->tIns_flnks[chang][GT]+=1;}
  else if(pre=='G' && post=='G') {T->tIns_flnks[chang][GG]+=1;}
  else if(pre=='G' && post=='C') {T->tIns_flnks[chang][GC]+=1;}
  else if(pre=='C' && post=='A') {T->tIns_flnks[chang][CA]+=1;}
  else if(pre=='C' && post=='T') {T->tIns_flnks[chang][CT]+=1;}
  else if(pre=='C' && post=='G') {T->tIns_flnks[chang][CG]+=1;}
  else if(pre=='C' && post=='C') {T->tIns_flnks[chang][CC]+=1;}
  }

}

/////////////////////new//////////////////////////////////
void get_seqInfo(bam1_t *b, char *ref, int ref_len, log_v *T, log_v *OA, int read_ovrlp_st, int read_ovrlp_nd)
 {
	 
         uint8_t *seq = bam1_seq(b), *qual = bam1_qual(b);
         int rf_pos=b->core.pos,rd_pos=0, seq_len=0, i,k,x,flg,pp;
         int qual_sum=0,sc_3p=0,sc_5p=0,hc_3p=0,hc_5p=0,Ins_c=0,Ins_f=0,mm_cnt=0;
	 int LOG[12], Acnt=0, Tcnt=0, Gcnt=0, Ccnt=0, Ncnt,iSize=0;
	 char CLP[200];
         uint32_t *cigar = bam_get_cigar(b);
         //////////////////Added on 25Sep202///////////////////////
         T->glog.tot_rd_cnt++;
         if( b->core.flag & BAM_FPAIRED) {T->glog.map_pr++;}
         if( b->core.flag & BAM_FPROPER_PAIR) {T->glog.propr_pr++;}
         //if( b->core.flag & BAM_FUNMAP) {T->glog.unmap++;}
         if( b->core.flag & BAM_FREVERSE) {T->glog.rev_cnt++;}
         if( b->core.flag & BAM_FREAD1) {T->glog.r1_cnt++;}
         if( b->core.flag & BAM_FREAD2) {T->glog.r2_cnt++;}
         if( b->core.flag & BAM_FSECONDARY) {T->glog.scndry_cnt++;}
         if( b->core.flag & BAM_FQCFAIL) {T->glog.qc_fail++;}
         if( b->core.flag & BAM_FSUPPLEMENTARY) {T->glog.supply_cnt++;}
         if((b->core.flag & BAM_FDUP)) {T->glog.DupCnt++;}
         if (b->core.tid != b->core.mtid) {T->glog.map_difchr++;}
         if((b->core.flag & BAM_FREVERSE) && (b->core.flag & BAM_FMREVERSE)) {T->glog.RR++;}
         else if(!(b->core.flag & BAM_FREVERSE) && !(b->core.flag & BAM_FMREVERSE)) {T->glog.FF++;}
         else {T->glog.FR++;}
         if(!(b->core.flag & BAM_FREVERSE) && (b->core.flag & BAM_FMREVERSE) && (b->core.mpos - b->core.pos) < -10 && (b->core.tid == b->core.mtid)) {T->glog.OO++;}
         else if((b->core.flag & BAM_FREVERSE) && !(b->core.flag & BAM_FMREVERSE) && (b->core.pos - b->core.mpos) < -10 && (b->core.tid == b->core.mtid)) {T->glog.OO++;}
            
         Acnt=0;Tcnt=0;Gcnt=0;Ccnt=0;
	 for(i=0;i<b->core.l_qseq;i++) 
	 {
               if(bam_nt16_rev_table[bam1_seqi(seq,i)]=='A') {Acnt+=1; if(!(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))) {T->glog.A_con[i]+=1;}}
          else if(bam_nt16_rev_table[bam1_seqi(seq,i)]=='T') {Tcnt+=1; if(!(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))) {T->glog.T_con[i]+=1;}}
          else if(bam_nt16_rev_table[bam1_seqi(seq,i)]=='G') {Gcnt+=1; if(!(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))) {T->glog.G_con[i]+=1;}}
          else if(bam_nt16_rev_table[bam1_seqi(seq,i)]=='C') {Ccnt+=1; if(!(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))) {T->glog.C_con[i]+=1;}}
	  
	  
          if(qual[i] >= MAX_BASEQ) {OA->glog.BaseQ_dist[MAX_BASEQ]+=1;} else {OA->glog.BaseQ_dist[qual[i]]+=1;}
	  
          
	 }
         
         if(!(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))) {OA->glog.GC_pct[(int)((float)(Gcnt+Ccnt)/(float)(Acnt+Tcnt+Gcnt+Ccnt)*100)]+=1;}
         
         
         for(i=0;i<b->core.l_qseq;i++) // No effect on bed overlap
	 {
           if(!(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))) { T->glog.Percycle_mn_baseQ[i]+=qual[i];}
         }
         
         
         
         if(b->core.isize < 0) {iSize=(b->core.isize * (-1));} else {iSize=b->core.isize;}
	 
	 if(iSize > MAX_ISIZE) {OA->glog.Isize_dist[MAX_ISIZE] += 1;} else {OA->glog.Isize_dist[iSize]+=1;}
	 
	 if(iSize >= 1000) {T->glog.Isize_ge1K++;}
	 
	 if(b->core.l_qseq > MAXREADLENGTH) {OA->glog.Rdlen[MAXREADLENGTH] += 1;} else {OA->glog.Rdlen[b->core.l_qseq]+=1;}
	 //iSize > MAX_ISIZE ? T->glog.Isize_dist[MAX_ISIZE]+=1; : T->glog.Isize_dist[MAX_ISIZE]+=1;
         
         if(b->core.qual >= MAX_MAPQ) {OA->glog.MapQ_dist[MAX_MAPQ] += 1;} else {OA->glog.MapQ_dist[b->core.qual]+=1;}
         
         if(b->core.qual == 0) {OA->glog.MapQ_bin_dist[1] = OA->glog.MapQ_bin_dist[1] + 1; }
         else if(b->core.qual <= 10) {OA->glog.MapQ_bin_dist[2] = OA->glog.MapQ_bin_dist[2] + 1;}
         else if(b->core.qual <= 20) {OA->glog.MapQ_bin_dist[3] = OA->glog.MapQ_bin_dist[3] + 1;}
         else if(b->core.qual <= 30) {OA->glog.MapQ_bin_dist[4] = OA->glog.MapQ_bin_dist[4] + 1;}
         else if(b->core.qual <= 40) {OA->glog.MapQ_bin_dist[5] = OA->glog.MapQ_bin_dist[5] + 1;}
         else if(b->core.qual <= 50) {OA->glog.MapQ_bin_dist[6] = OA->glog.MapQ_bin_dist[6] + 1;}
         else {OA->glog.MapQ_bin_dist[7] = OA->glog.MapQ_bin_dist[7] + 1;}
         //////////////////////////////////////////////////////////
  
         for(k=0;k<12;k++){LOG[k]=0;}
    
         for (k = 0 ; k < b->core.n_cigar; ++k) 
         {
             int l = cigar[k]>>4;
             int op = cigar[k]&0xf;
             if(op == BAM_CDEL) {T->midc.del_cnt++; T->midc.del_cnt_bs+=l;}
             if(op == BAM_CINS) {T->midc.ins_cnt++; T->midc.ins_cnt_bs+=l;}
             if(op == BAM_CHARD_CLIP) {T->midc.hrdclp++; T->midc.hrdclp_bs+=l;}
             if(op == BAM_CSOFT_CLIP) 
             {
                 T->midc.sftclp++; T->midc.sftclp_bs+=l;
                 x=0;CLP[0]='\0';
                 if(l>=12)
                 {
                     for(i=rd_pos;i<rd_pos+l;i++) {CLP[x]=bam_nt16_rev_table[bam1_seqi(seq,i)];x++;}
                     CLP[x]='\0';
                     x=0;
                     while(x < 6 )
                     {
                          for(i=0;i<l-11;i++)
                          {
                               flg=1;
                               for(pp=i;pp<(i+11) && flg;pp++)
                               {
                                    if(CLP[pp] != Adp_prim_trans[x][pp-i]) {flg=0;}
                               }
                              
                               if(flg==1)
                               {
                                    if(x==0) {T->adptr.NTS++;}
                                    else if(x==1) {T->adptr.IPCRprm++;}
                                    else {T->adptr.IADptr++;}
                               }
                          }
                          x++;
                     }
                 }
             }
             
             if(op == BAM_CHARD_CLIP && !k) 
             {
                 hc_5p=l; 
                 if(l>20) { T->tClp3[4]+=1;}
                 else if(l>=16 && l<=20) { T->tClp3[3]+=1; }
                 else if(l>=11 && l<=15) { T->tClp3[2]+=1; }
                 else if(l>= 6 && l<=10) { T->tClp3[1]+=1; }
                 else { T->tClp3[0]+=1; }
                 seq_len+=l;
             }
             else if(op == BAM_CHARD_CLIP) 
             {
                 hc_3p=l; 
                 if(l>20) { T->tClp5[4]+=1; }
                 else if(l>=16 && l<=20) { T->tClp5[3]+=1; }
                 else if(l>=11 && l<=15) { T->tClp5[2]+=1; }
                 else if(l>= 6 && l<=10) { T->tClp5[1]+=1; }
                 else { T->tClp5[0]+=1; } 
                 seq_len+=l;
             }
             else if(op == BAM_CSOFT_CLIP && !k)
             {
                 sc_5p=l;
                 if(rd_pos >= read_ovrlp_st && rd_pos <= read_ovrlp_nd)
                 {
                 if(l>20) { T->tClp3[4]+=1; }
                 else if(l>=16 && l<=20) { T->tClp3[3]+=1; }
                 else if(l>=11 && l<=15) { T->tClp3[2]+=1; }
                 else if(l>= 6 && l<=10) { T->tClp3[1]+=1; }
                 else { T->tClp3[0]+=1; }
                 }
                 rd_pos+=l;
             }
             else if(op == BAM_CSOFT_CLIP) 
             {
                 sc_3p=l;
                 if(rd_pos >= read_ovrlp_st && rd_pos <= read_ovrlp_nd)
                 {
                 if(l>20) { T->tClp5[4]+=1; }
                 else if(l>=16 && l<=20) { T->tClp5[3]+=1; }
                 else if(l>=11 && l<=15) { T->tClp5[2]+=1; }
                 else if(l>= 6 && l<=10) { T->tClp5[1]+=1; }
                 else { T->tClp5[0]+=1; }
                 }
                 rd_pos+=l;
             }
             else if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CBACK)
             {
                for(i=0;i<l;i++)
                {
                  if(rd_pos+i >= read_ovrlp_st && rd_pos+i <= read_ovrlp_nd)
                  {
                    if(bam_nt16_rev_table[bam1_seqi(seq, rd_pos+i)]=='A') { T->glog.AAcnt+=1; }
                    else if(bam_nt16_rev_table[bam1_seqi(seq, rd_pos+i)]=='T') { T->glog.TTcnt+=1; }
                    else if(bam_nt16_rev_table[bam1_seqi(seq, rd_pos+i)]=='G') { T->glog.GGcnt+=1; }
                    else if(bam_nt16_rev_table[bam1_seqi(seq, rd_pos+i)]=='C') { T->glog.CCcnt+=1; }
	            else { T->glog.NNcnt++; }
                    
                   if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)] != ((ref && (rf_pos+i) < ref_len)? toupper(ref[rf_pos+i]) : 'N'))
                   {
//printf("Base=%c	ref=%c\n",bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)], toupper(ref[rf_pos+i]));
                        if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]!='N' && toupper(ref[rf_pos+i])!='N') {mm_cnt+=1;
                        //if(b->core.flag & BAM_FREAD1) { if(b->core.flag & BAM_FREVERSE) {printf("Read1_rs	");} else{printf("Read1_fs	");}}
           		 //else { if(b->core.flag & BAM_FREVERSE) {printf("Read2_rs	");} else{printf("Read2_fs	");}}
                        //printf("MMBC_Q_context	%c>%c	%c	%c	%d	%d\n", toupper(ref[rf_pos+i]), bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)], toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), qual[rd_pos+i], rd_pos+i);
                        }
                        if(toupper(ref[rf_pos+i])=='A' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='T')      
                        {T->MMC[AT]+=1; count_flanks(toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), AT, 0, T); binQual(qual[rd_pos+i], AT, T);}
                        else if(toupper(ref[rf_pos+i])=='A' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='C') 
                        {T->MMC[AC]+=1; count_flanks(toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), AC, 0, T); binQual(qual[rd_pos+i], AC, T);}
                        else if(toupper(ref[rf_pos+i])=='A' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='G') 
                        {T->MMC[AG]+=1; count_flanks(toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), AG, 0, T); binQual(qual[rd_pos+i], AG, T);}
                        else if(toupper(ref[rf_pos+i])=='T' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='A') 
                        {T->MMC[TA]+=1; count_flanks(toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), TA, 0, T); binQual(qual[rd_pos+i], TA, T);}
                        else if(toupper(ref[rf_pos+i])=='T' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='C') 
                        {T->MMC[TC]+=1; count_flanks(toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), TC, 0, T); binQual(qual[rd_pos+i], TC, T);}
                        else if(toupper(ref[rf_pos+i])=='T' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='G') 
                        {T->MMC[TG]+=1; count_flanks(toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), TG, 0, T); binQual(qual[rd_pos+i], TG, T);}
                        else if(toupper(ref[rf_pos+i])=='C' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='A') 
                        {T->MMC[CA]+=1; count_flanks(toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), CA, 0, T); binQual(qual[rd_pos+i], CA, T);}
                        else if(toupper(ref[rf_pos+i])=='C' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='T') 
                        {T->MMC[CT]+=1; count_flanks(toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), CT, 0, T); binQual(qual[rd_pos+i], CT, T);}
                        else if(toupper(ref[rf_pos+i])=='C' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='G') 
                        {T->MMC[CG]+=1; count_flanks(toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), CG, 0, T); binQual(qual[rd_pos+i], CG, T);}
                        else if(toupper(ref[rf_pos+i])=='G' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='A') 
                        {T->MMC[GA]+=1; count_flanks(toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), GA, 0, T); binQual(qual[rd_pos+i], GA, T);}
                        else if(toupper(ref[rf_pos+i])=='G' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='T') 
                        {T->MMC[GT]+=1; count_flanks(toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), GT, 0, T); binQual(qual[rd_pos+i], GT, T);}
                        else if(toupper(ref[rf_pos+i])=='G' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='C') 
                        {T->MMC[GC]+=1; count_flanks(toupper(ref[(rf_pos+i)-1]), toupper(ref[(rf_pos+i)+1]), GC, 0, T); binQual(qual[rd_pos+i], GC, T);}     
                   }
                 }
                }
                rf_pos+=l; rd_pos+=l;
             }
             else if(op == BAM_CINS) 
             {
                Ins_c+=l;
                if(l==1 && (rd_pos >= read_ovrlp_st && rd_pos <= read_ovrlp_nd))
                {
                   if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='A' && !(T->Ins1.A)) {T->Ins1.A+=1;}
                   else if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='T' && !(T->Ins1.T)) {T->Ins1.T+=1;}
                   else if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='G' && !(T->Ins1.G)) {T->Ins1.G+=1;}
                   else if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='C' && !(T->Ins1.C)) {T->Ins1.C+=1;}
                   count_flanks(toupper(ref[(rf_pos-1)]), toupper(ref[rf_pos]), 0, 2, T);
                   T->tIns[0]+=1;
                        if((bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='A' && toupper(ref[rf_pos-2])=='A' && toupper(ref[rf_pos-1])=='A' && toupper(ref[rf_pos])=='A') || (bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='A' && toupper(ref[rf_pos+2])=='A' && toupper(ref[rf_pos+1])=='A' && toupper(ref[rf_pos])=='A')) {T->midc.polyIns++;}
                   else if((bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='T' && toupper(ref[rf_pos-2])=='T' && toupper(ref[rf_pos-1])=='T' && toupper(ref[rf_pos])=='T') || (bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='T' && toupper(ref[rf_pos+2])=='T' && toupper(ref[rf_pos+1])=='T' && toupper(ref[rf_pos])=='T')) {T->midc.polyIns++;}
                   else if((bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='G' && toupper(ref[rf_pos-2])=='G' && toupper(ref[rf_pos-1])=='G' && toupper(ref[rf_pos])=='G') || (bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='G' && toupper(ref[rf_pos+2])=='G' && toupper(ref[rf_pos+1])=='G' && toupper(ref[rf_pos])=='G')) {T->midc.polyIns++;}
                   else if((bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='C' && toupper(ref[rf_pos-2])=='C' && toupper(ref[rf_pos-1])=='C' && toupper(ref[rf_pos])=='C') || (bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='C' && toupper(ref[rf_pos+2])=='C' && toupper(ref[rf_pos+1])=='C' && toupper(ref[rf_pos])=='C')) {T->midc.polyIns++;}
                }
                else if(l>=5)  {if(!(T->tIns[4]) && (rd_pos >= read_ovrlp_st && rd_pos <= read_ovrlp_nd)) {T->tIns[4]+=1; count_flanks(toupper(ref[(rf_pos-1)]), toupper(ref[rf_pos]), 4, 2, T);}} 
                else {if(!(T->tIns[l-1]) && (rd_pos >= read_ovrlp_st && rd_pos <= read_ovrlp_nd)) {T->tIns[l-1]+=1; count_flanks(toupper(ref[(rf_pos-1)]), toupper(ref[rf_pos]), (l-1), 2, T);}}
                
                rd_pos+=l;
             }
             else if(op == BAM_CDEL)
             {
                if(l==1 && (rd_pos >= read_ovrlp_st && rd_pos <= read_ovrlp_nd))
                {
                   // The base should be get from reference
                   if(toupper(ref[rf_pos])=='A' && !(T->Del1.A)) {T->Del1.A+=1;}
                   else if(toupper(ref[rf_pos])=='T' && !(T->Del1.T)) {T->Del1.T+=1;}
                   else if(toupper(ref[rf_pos])=='G' && !(T->Del1.G)) {T->Del1.G+=1;}
                   else if(toupper(ref[rf_pos])=='C' && !(T->Del1.C)) {T->Del1.C+=1;}
                   count_flanks(toupper(ref[rf_pos-1]), toupper(ref[(rf_pos+l)]), 0, 1, T);
                   if((toupper(ref[rf_pos])=='A' && toupper(ref[rf_pos-3])=='A' && toupper(ref[rf_pos-2])=='A' && toupper(ref[rf_pos-1])=='A') || (toupper(ref[rf_pos])=='A' && toupper(ref[rf_pos+3])=='A' && toupper(ref[rf_pos+2])=='A' && toupper(ref[rf_pos+1])=='A')) {T->midc.polyDel++;}
                   else if((toupper(ref[rf_pos])=='T' && toupper(ref[rf_pos-3])=='T' && toupper(ref[rf_pos-2])=='T' && toupper(ref[rf_pos-1])=='T') || (toupper(ref[rf_pos])=='T' && toupper(ref[rf_pos+3])=='T' && toupper(ref[rf_pos+2])=='T' && toupper(ref[rf_pos+1])=='T')) {T->midc.polyDel++;}
                   else if((toupper(ref[rf_pos])=='G' && toupper(ref[rf_pos-3])=='G' && toupper(ref[rf_pos-2])=='G' && toupper(ref[rf_pos-1])=='G') || (toupper(ref[rf_pos])=='G' && toupper(ref[rf_pos+3])=='G' && toupper(ref[rf_pos+2])=='G' && toupper(ref[rf_pos+1])=='G')) {T->midc.polyDel++;}
                   else if((toupper(ref[rf_pos])=='C' && toupper(ref[rf_pos-3])=='C' && toupper(ref[rf_pos-2])=='C' && toupper(ref[rf_pos-1])=='C') || (toupper(ref[rf_pos])=='C' && toupper(ref[rf_pos+3])=='C' && toupper(ref[rf_pos+2])=='C' && toupper(ref[rf_pos+1])=='C')) {T->midc.polyDel++;}
                   T->tDel[0]+=1;
                   /*if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='A' && !(T->Del1.A)) {T->Del1.A+=1;}
                   else if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='T' && !(T->Del1.T)) {T->Del1.T+=1;}
                   else if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='G' && !(T->Del1.G)) {T->Del1.G+=1;}
                   else if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+1)]=='C' && !(T->Del1.C)) {T->Del1.C+=1;}*/
                }
                else if(l>=5) {if(!(T->tDel[4]) && (rd_pos >= read_ovrlp_st && rd_pos <= read_ovrlp_nd)) {T->tDel[4]+=1; count_flanks(toupper(ref[rf_pos-1]), toupper(ref[(rf_pos+l)]), 4, 1, T);}} 
                else {if(!(T->tDel[l-1]) && (rd_pos >= read_ovrlp_st && rd_pos <= read_ovrlp_nd)) {T->tDel[l-1]+=1; count_flanks(toupper(ref[rf_pos-1]), toupper(ref[(rf_pos+l)]), (l-1), 1, T);}}
                rf_pos+=l;
             }
             else if(op == BAM_CREF_SKIP) {rf_pos+=l;}
             
         }

	 //printf("MMC = %d	MD = %s\n",mm_cnt,bam_aux_get(b, "MD"));
         //if(mm_cnt>=5) {if(!(T->tmmc[4])) {T->tmmc[4]+=1;}} else if(mm_cnt>=1) {if(!(T->tmmc[mm_cnt-1])) {T->tmmc[mm_cnt-1]+=1;}}
         
         if(mm_cnt) { T->midc.basemmc+=mm_cnt; if(mm_cnt>=5) {T->tmmc[4]+=1;} else {T->tmmc[mm_cnt-1]+=1;}}
         Ins_f=hc_3p+hc_5p+sc_3p+sc_5p+Ins_c;
         //printf("HC3=%d	HC5=%d	SC3=%d	SC5=%d	Insertion=%d\n",hc_3p,hc_5p,sc_3p,sc_5p,Ins_c);
         rd_pos=0;seq_len=0;rf_pos=b->core.pos; 
         //if(Ins_f==0) {T->no_indclp++;}
         for (k = 0 ; k < b->core.n_cigar && Ins_f==0 ; ++k) 
         {
                
                int l = cigar[k]>>4;
                int op = cigar[k]&0xf;
                if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF )
                {
                   for(i=0;i<l;i++)
                   {
                       if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)] != ((ref && (rf_pos+i) < ref_len)? toupper(ref[rf_pos+i]) : 'N'))
                       {
                               //printf("%c>%c[Readpos=%d]	BaseQual=%d\n",ref[rf_pos+i],bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)],rd_pos+i,qual[rd_pos+i]);
                               if(!T->no_indclp) {T->no_indclp++;};
                               T->rmmQV.rd_posi[rd_pos+i]+=qual[rd_pos+i];
                               if(toupper(ref[rf_pos+i])=='A' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='T')      {T->rmmc[AT].rd_posi[rd_pos+i]+=1; /*T->rmmQV[AT]+=qual[rd_pos+i];*/}
                               else if(toupper(ref[rf_pos+i])=='A' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='C') {T->rmmc[AC].rd_posi[rd_pos+i]+=1; /*T->rmmQV[AC]+=qual[rd_pos+i];*/}
                               else if(toupper(ref[rf_pos+i])=='A' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='G') {T->rmmc[AG].rd_posi[rd_pos+i]+=1; /*T->rmmQV[AG]+=qual[rd_pos+i];*/}
                               else if(toupper(ref[rf_pos+i])=='T' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='A') {T->rmmc[TA].rd_posi[rd_pos+i]+=1; /*T->rmmQV[TA]+=qual[rd_pos+i];*/}
                               else if(toupper(ref[rf_pos+i])=='T' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='C') {T->rmmc[TC].rd_posi[rd_pos+i]+=1; /*T->rmmQV[TC]+=qual[rd_pos+i];*/}
                               else if(toupper(ref[rf_pos+i])=='T' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='G') {T->rmmc[TG].rd_posi[rd_pos+i]+=1; /*T->rmmQV[TG]+=qual[rd_pos+i];*/}
                               else if(toupper(ref[rf_pos+i])=='C' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='A') {T->rmmc[CA].rd_posi[rd_pos+i]+=1; /*T->rmmQV[CA]+=qual[rd_pos+i];*/}
                               else if(toupper(ref[rf_pos+i])=='C' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='T') {T->rmmc[CT].rd_posi[rd_pos+i]+=1; /*T->rmmQV[CT]+=qual[rd_pos+i];*/}
                               else if(toupper(ref[rf_pos+i])=='C' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='G') {T->rmmc[CG].rd_posi[rd_pos+i]+=1; /*T->rmmQV[CG]+=qual[rd_pos+i];*/}
                               else if(toupper(ref[rf_pos+i])=='G' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='A') {T->rmmc[GA].rd_posi[rd_pos+i]+=1; /*T->rmmQV[GA]+=qual[rd_pos+i];*/}
                               else if(toupper(ref[rf_pos+i])=='G' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='T') {T->rmmc[GT].rd_posi[rd_pos+i]+=1; /*T->rmmQV[GT]+=qual[rd_pos+i];*/}
                               else if(toupper(ref[rf_pos+i])=='G' && bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]=='C') {T->rmmc[GC].rd_posi[rd_pos+i]+=1; /*T->rmmQV[GC]+=qual[rd_pos+i];*/}
                       }
                   }
                   rf_pos+=l; rd_pos+=l;
                }
                else if(op == BAM_CDEL || op == BAM_CREF_SKIP)
                {rf_pos+=l;}
                
         }
            
         
 }

void PrintLogFiles(log_v *OA, int type, int rLen, char *outdir)
{
      FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8, *fp9, *fp10;
      FILE *fp1c;
      char *outfn;
      outfn = calloc(strlen(outdir) + 500, 1);
      if(type == 1)
      {
         *outfn=0;
         sprintf(outfn, "%s/Overall_SN_MM.txt", outdir);
         fp1=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Overall_SN_MM_count.txt", outdir);
         fp1c=fopen(outfn,"w");
         
         *outfn=0;
         sprintf(outfn, "%s/Overall_INDEL_MM.txt", outdir);
         fp2=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Overall_CLIP_MM.txt", outdir);
         fp3=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Overall_BaseChange_SN_MM.txt", outdir);
         fp4=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Overall_BaseChangeQV_SN_MM.txt", outdir);
         fp5=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Overall_ReadposBaseChange_SN_MM.txt", outdir);
         fp6=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Overall_ReadposTotal_SN_MM.txt", outdir);
         fp7=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Overall_Readpos_MMQV.txt", outdir);
         fp8=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Overall_BaseChange_SN_MM_flanks.txt", outdir);
         fp9=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Overall_INDEL_MM_flanks.txt", outdir);
         fp10=fopen(outfn,"w");
      }
      if(type == 2)
      {
         *outfn=0;
         sprintf(outfn, "%s/Read1_fs_SN_MM.txt", outdir);
         fp1=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_fs_SN_MM_count.txt", outdir);
         fp1c=fopen(outfn,"w");
         
         *outfn=0;
         sprintf(outfn, "%s/Read1_fs_INDEL_MM.txt", outdir);
         fp2=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_fs_CLIP_MM.txt", outdir);
         fp3=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_fs_BaseChange_SN_MM.txt", outdir);
         fp4=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_fs_BaseChangeQV_SN_MM.txt", outdir);
         fp5=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_fs_ReadposBaseChange_SN_MM.txt", outdir);
         fp6=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_fs_ReadposTotal_SN_MM.txt", outdir);
         fp7=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_fs_Readpos_MMQV.txt", outdir);
         fp8=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_fs_BaseChange_SN_MM_flanks.txt", outdir);
         fp9=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_fs_INDEL_MM_flanks.txt", outdir);
         fp10=fopen(outfn,"w");
      }
      if(type == 3)
      {
         *outfn=0;
         sprintf(outfn, "%s/Read1_rs_SN_MM.txt", outdir);
         fp1=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_rs_SN_MM_count.txt", outdir);
         fp1c=fopen(outfn,"w");
         
         *outfn=0;
         sprintf(outfn, "%s/Read1_rs_INDEL_MM.txt", outdir);
         fp2=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_rs_CLIP_MM.txt", outdir);
         fp3=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_rs_BaseChange_SN_MM.txt", outdir);
         fp4=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_rs_BaseChangeQV_SN_MM.txt", outdir);
         fp5=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_rs_ReadposBaseChange_SN_MM.txt", outdir);
         fp6=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_rs_ReadposTotal_SN_MM.txt", outdir);
         fp7=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_rs_Readpos_MMQV.txt", outdir);
         fp8=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_rs_BaseChange_SN_MM_flanks.txt", outdir);
         fp9=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read1_rs_INDEL_MM_flanks.txt", outdir);
         fp10=fopen(outfn,"w");
      }
      if(type == 4)
      {
         *outfn=0;
         sprintf(outfn, "%s/Read2_fs_SN_MM.txt", outdir);
         fp1=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_fs_SN_MM_count.txt", outdir);
         fp1c=fopen(outfn,"w");
         
         *outfn=0;
         sprintf(outfn, "%s/Read2_fs_INDEL_MM.txt", outdir);
         fp2=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_fs_CLIP_MM.txt", outdir);
         fp3=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_fs_BaseChange_SN_MM.txt", outdir);
         fp4=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_fs_BaseChangeQV_SN_MM.txt", outdir);
         fp5=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_fs_ReadposBaseChange_SN_MM.txt", outdir);
         fp6=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_fs_ReadposTotal_SN_MM.txt", outdir);
         fp7=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_fs_Readpos_MMQV.txt", outdir);
         fp8=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_fs_BaseChange_SN_MM_flanks.txt", outdir);
         fp9=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_fs_INDEL_MM_flanks.txt", outdir);
         fp10=fopen(outfn,"w");
      }
      if(type == 5)
      {
         *outfn=0;
         sprintf(outfn, "%s/Read2_rs_SN_MM.txt", outdir);
         fp1=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_rs_SN_MM_count.txt", outdir);
         fp1c=fopen(outfn,"w");
         
         *outfn=0;
         sprintf(outfn, "%s/Read2_rs_INDEL_MM.txt", outdir);
         fp2=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_rs_CLIP_MM.txt", outdir);
         fp3=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_rs_BaseChange_SN_MM.txt", outdir);
         fp4=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_rs_BaseChangeQV_SN_MM.txt", outdir);
         fp5=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_rs_ReadposBaseChange_SN_MM.txt", outdir);
         fp6=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_rs_ReadposTotal_SN_MM.txt", outdir);
         fp7=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_rs_Readpos_MMQV.txt", outdir);
         fp8=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_rs_BaseChange_SN_MM_flanks.txt", outdir);
         fp9=fopen(outfn,"w");
         *outfn=0;
         sprintf(outfn, "%s/Read2_rs_INDEL_MM_flanks.txt", outdir);
         fp10=fopen(outfn,"w");
      }
      int kk, ii, rdpos_total[rLen];
      long int totl=0;
      //printf("Miss-matchs logs::\n");
      
      fprintf(fp1,"Missmatch	pct_of_reads\n");
      for(kk=0;kk<4;kk++) 
      {fprintf(fp1,"%d	%lf\n",kk+1,((OA->tmmc[kk]/(double)OA->cnt)*100));}
      fprintf(fp1,">=5	%lf\n",((OA->tmmc[4]/(double)OA->cnt)*100));
      
      ////*updated on 19Apr2021 start*///
      fprintf(fp1c,"Missmatch	no_of_reads\n");
      fprintf(fp1c,"0	%ld\n",(OA->cnt - (OA->tmmc[0]+OA->tmmc[1]+OA->tmmc[2]+OA->tmmc[3]+OA->tmmc[4])));
      for(kk=0;kk<4;kk++)
      {fprintf(fp1c,"%d	%ld\n",kk+1, OA->tmmc[kk]);}
      fprintf(fp1c,">=5	%ld\n",OA->tmmc[4]);
      ////*updated on 19Apr2021 end*///
      //printf("\n");
      fprintf(fp2,"Len	Type	pct_of_reads\n");
      fprintf(fp2,"0	A	0\n");
      fprintf(fp2,"-1	A	%lf\n",((OA->Del1.A/(double)OA->cnt)*100));
      fprintf(fp2,"-1	T	%lf\n",((OA->Del1.T/(double)OA->cnt)*100));
      fprintf(fp2,"-1	G	%lf\n",((OA->Del1.G/(double)OA->cnt)*100));
      fprintf(fp2,"-1	C	%lf\n",((OA->Del1.C/(double)OA->cnt)*100));
      for(kk=1;kk<4;kk++)
      {
         fprintf(fp2,"-%d	Multi-nucleotide	%lf\n",kk+1,((OA->tDel[kk]/(double)OA->cnt)*100));
         fprintf(fp2,"-%d	A	0.0\n",kk+1);
         fprintf(fp2,"-%d	T	0.0\n",kk+1);
         fprintf(fp2,"-%d	G	0.0\n",kk+1);
         fprintf(fp2,"-%d	C	0.0\n",kk+1);
      }
      fprintf(fp2,">=-5	Multi-nucleotide	%lf\n",((OA->tDel[4]/(double)OA->cnt)*100));
      fprintf(fp2,">=-5	A	0.0\n");
      fprintf(fp2,">=-5	T	0.0\n");
      fprintf(fp2,">=-5	G	0.0\n");
      fprintf(fp2,">=-5	C	0.0\n");

      fprintf(fp2,"+1	A	%f\n",((OA->Ins1.A/(double)OA->cnt)*100));
      fprintf(fp2,"+1	T	%f\n",((OA->Ins1.T/(double)OA->cnt)*100));
      fprintf(fp2,"+1	G	%f\n",((OA->Ins1.G/(double)OA->cnt)*100));
      fprintf(fp2,"+1	C	%f\n",((OA->Ins1.C/(double)OA->cnt)*100));
      for(kk=1;kk<4;kk++)
      {
         fprintf(fp2,"+%d	Multi-nucleotide	%f\n",kk+1,((OA->tIns[kk]/(double)OA->cnt)*100));
         fprintf(fp2,"+%d	A	0.0\n",kk+1);
         fprintf(fp2,"+%d	T	0.0\n",kk+1);
         fprintf(fp2,"+%d	G	0.0\n",kk+1);
         fprintf(fp2,"+%d	C	0.0\n",kk+1);
      }
      fprintf(fp2,">=+5	Multi-nucleotide	%f\n",((OA->tIns[4]/(double)OA->cnt)*100));
      fprintf(fp2,">=+5	A	0.0\n");
      fprintf(fp2,">=+5	T	0.0\n");
      fprintf(fp2,">=+5	G	0.0\n");
      fprintf(fp2,">=+5	C	0.0\n");
      
      
      fprintf(fp3,"Bases	End	pct_of_reads\n");
      fprintf(fp3,"1-5	5prime	%f\n",((OA->tClp5[0]/(double)OA->cnt)*100));
      fprintf(fp3,"6-10	5prime	%f\n",((OA->tClp5[1]/(double)OA->cnt)*100));
      fprintf(fp3,"11-15	5prime	%f\n",((OA->tClp5[2]/(double)OA->cnt)*100));
      fprintf(fp3,"16-20	5prime	%f\n",((OA->tClp5[3]/(double)OA->cnt)*100));
      fprintf(fp3,">20	5prime	%f\n",((OA->tClp5[4]/(double)OA->cnt)*100));
      fprintf(fp3,"1-5	3prime	%f\n",((OA->tClp3[0]/(double)OA->cnt)*100));
      fprintf(fp3,"6-10	3prime	%f\n",((OA->tClp3[1]/(double)OA->cnt)*100));
      fprintf(fp3,"11-15	3prime	%f\n",((OA->tClp3[2]/(double)OA->cnt)*100));
      fprintf(fp3,"16-20	3prime	%f\n",((OA->tClp3[3]/(double)OA->cnt)*100));
      fprintf(fp3,">20	3prime	%f\n",((OA->tClp3[4]/(double)OA->cnt)*100));
      
      for(kk=0;kk<12;kk++) 
      {
         totl+=OA->MMC[kk];
      }
      //printf("\n");
      fprintf(fp4,"Change	pct_of_bases\n");

      for(kk=0;kk<12;kk++) 
      {
         if(kk==AT) {fprintf(fp4,"A>T	%f\n",((OA->MMC[AT]/(float)totl)*100));}
         else if(kk==AC) {fprintf(fp4,"A>C	%f\n",((OA->MMC[AC]/(float)totl)*100));}
         else if(kk==AG) {fprintf(fp4,"A>G	%f\n",((OA->MMC[AG]/(float)totl)*100));}
         else if(kk==TA) {fprintf(fp4,"T>A	%f\n",((OA->MMC[TA]/(float)totl)*100));}
         else if(kk==TC) {fprintf(fp4,"T>C	%f\n",((OA->MMC[TC]/(float)totl)*100));}
         else if(kk==TG) {fprintf(fp4,"T>G	%f\n",((OA->MMC[TG]/(float)totl)*100));}
         else if(kk==CA) {fprintf(fp4,"C>A	%f\n",((OA->MMC[CA]/(float)totl)*100));}
         else if(kk==CT) {fprintf(fp4,"C>T	%f\n",((OA->MMC[CT]/(float)totl)*100));}
         else if(kk==CG) {fprintf(fp4,"C>G	%f\n",((OA->MMC[CG]/(float)totl)*100));}
         else if(kk==GA) {fprintf(fp4,"G>A	%f\n",((OA->MMC[GA]/(float)totl)*100));}
         else if(kk==GT) {fprintf(fp4,"G>T	%f\n",((OA->MMC[GT]/(float)totl)*100));}
         else if(kk==GC) {fprintf(fp4,"G>C	%f\n",((OA->MMC[GC]/(float)totl)*100));}
      }
      //printf("\n");
      fprintf(fp5,"Change	Range	QV\n");
      for(kk=0;kk<12;kk++)
      {
        totl=OA->mmqv[kk].qv1+OA->mmqv[kk].qv2+OA->mmqv[kk].qv3+OA->mmqv[kk].qv8;
         if(kk==AT) {fprintf(fp5,"A>T	<=10	%f\nA>T	11-20	%f\nA>T	21-30	%f\nA>T	>30	%f\n",((OA->mmqv[kk].qv1/(float)totl)*100),((OA->mmqv[kk].qv2/(float)totl)*100),((OA->mmqv[kk].qv3/(float)totl)*100),((OA->mmqv[kk].qv8/(float)totl)*100));}
         else if(kk==AC) {fprintf(fp5,"A>C	<=10	%f\nA>C	11-20	%f\nA>C	21-30	%f\nA>C	>30	%f\n",((OA->mmqv[kk].qv1/(float)totl)*100),((OA->mmqv[kk].qv2/(float)totl)*100),((OA->mmqv[kk].qv3/(float)totl)*100),((OA->mmqv[kk].qv8/(float)totl)*100));}
         else if(kk==AG) {fprintf(fp5,"A>G	<=10	%f\nA>G	11-20	%f\nA>G	21-30	%f\nA>G	>30	%f\n",((OA->mmqv[kk].qv1/(float)totl)*100),((OA->mmqv[kk].qv2/(float)totl)*100),((OA->mmqv[kk].qv3/(float)totl)*100),((OA->mmqv[kk].qv8/(float)totl)*100));}
         else if(kk==TA) {fprintf(fp5,"T>A	<=10	%f\nT>A	11-20	%f\nT>A	21-30	%f\nT>A	>30	%f\n",((OA->mmqv[kk].qv1/(float)totl)*100),((OA->mmqv[kk].qv2/(float)totl)*100),((OA->mmqv[kk].qv3/(float)totl)*100),((OA->mmqv[kk].qv8/(float)totl)*100));}
         else if(kk==TC) {fprintf(fp5,"T>C	<=10	%f\nT>C	11-20	%f\nT>C	21-30	%f\nT>C	>30	%f\n",((OA->mmqv[kk].qv1/(float)totl)*100),((OA->mmqv[kk].qv2/(float)totl)*100),((OA->mmqv[kk].qv3/(float)totl)*100),((OA->mmqv[kk].qv8/(float)totl)*100));}
         else if(kk==TG) {fprintf(fp5,"T>G	<=10	%f\nT>G	11-20	%f\nT>G	21-30	%f\nT>G	>30	%f\n",((OA->mmqv[kk].qv1/(float)totl)*100),((OA->mmqv[kk].qv2/(float)totl)*100),((OA->mmqv[kk].qv3/(float)totl)*100),((OA->mmqv[kk].qv8/(float)totl)*100));}
         else if(kk==CA) {fprintf(fp5,"C>A	<=10	%f\nC>A	11-20	%f\nC>A	21-30	%f\nC>A	>30	%f\n",((OA->mmqv[kk].qv1/(float)totl)*100),((OA->mmqv[kk].qv2/(float)totl)*100),((OA->mmqv[kk].qv3/(float)totl)*100),((OA->mmqv[kk].qv8/(float)totl)*100));}
         else if(kk==CT) {fprintf(fp5,"C>T	<=10	%f\nC>T	11-20	%f\nC>T	21-30	%f\nC>T	>30	%f\n",((OA->mmqv[kk].qv1/(float)totl)*100),((OA->mmqv[kk].qv2/(float)totl)*100),((OA->mmqv[kk].qv3/(float)totl)*100),((OA->mmqv[kk].qv8/(float)totl)*100));}
         else if(kk==CG) {fprintf(fp5,"C>G	<=10	%f\nC>G	11-20	%f\nC>G	21-30	%f\nC>G	>30	%f\n",((OA->mmqv[kk].qv1/(float)totl)*100),((OA->mmqv[kk].qv2/(float)totl)*100),((OA->mmqv[kk].qv3/(float)totl)*100),((OA->mmqv[kk].qv8/(float)totl)*100));}
         else if(kk==GA) {fprintf(fp5,"G>A	<=10	%f\nG>A	11-20	%f\nG>A	21-30	%f\nG>A	>30	%f\n",((OA->mmqv[kk].qv1/(float)totl)*100),((OA->mmqv[kk].qv2/(float)totl)*100),((OA->mmqv[kk].qv3/(float)totl)*100),((OA->mmqv[kk].qv8/(float)totl)*100));}
         else if(kk==GT) {fprintf(fp5,"G>T	<=10	%f\nG>T	11-20	%f\nG>T	21-30	%f\nG>T	>30	%f\n",((OA->mmqv[kk].qv1/(float)totl)*100),((OA->mmqv[kk].qv2/(float)totl)*100),((OA->mmqv[kk].qv3/(float)totl)*100),((OA->mmqv[kk].qv8/(float)totl)*100));}
         else if(kk==GC) {fprintf(fp5,"G>C	<=10	%f\nG>C	11-20	%f\nG>C	21-30	%f\nG>C	>30	%f\n",((OA->mmqv[kk].qv1/(float)totl)*100),((OA->mmqv[kk].qv2/(float)totl)*100),((OA->mmqv[kk].qv3/(float)totl)*100),((OA->mmqv[kk].qv8/(float)totl)*100));}
      }
  
      for(ii=0;ii<rLen;ii++)
      {
         rdpos_total[ii]=0;
         for(kk=0;kk<12;kk++)
         {
           rdpos_total[ii]+=OA->rmmc[kk].rd_posi[ii];
         }
      }
      //printf("\n");
      
      fprintf(fp6,"Change	Pos	Value\n");
      for(ii=0;ii<rLen;ii++)
      {
         for(kk=0;kk<12;kk++)
         {
            if(kk==AT) {fprintf(fp6,"A>T	%d	%f\n",ii+1,((OA->rmmc[AT].rd_posi[ii]/(float)rdpos_total[ii])*100));}
            else if(kk==AC) {fprintf(fp6,"A>C	%d	%f\n",ii+1,((OA->rmmc[AC].rd_posi[ii]/(float)rdpos_total[ii])*100));}
            else if(kk==AG) {fprintf(fp6,"A>G	%d	%f\n",ii+1,((OA->rmmc[AG].rd_posi[ii]/(float)rdpos_total[ii])*100));}
            else if(kk==TA) {fprintf(fp6,"T>A	%d	%f\n",ii+1,((OA->rmmc[TA].rd_posi[ii]/(float)rdpos_total[ii])*100));}
            else if(kk==TC) {fprintf(fp6,"T>C	%d	%f\n",ii+1,((OA->rmmc[TC].rd_posi[ii]/(float)rdpos_total[ii])*100));}
            else if(kk==TG) {fprintf(fp6,"T>G	%d	%f\n",ii+1,((OA->rmmc[TG].rd_posi[ii]/(float)rdpos_total[ii])*100));}
            else if(kk==CA) {fprintf(fp6,"C>A	%d	%f\n",ii+1,((OA->rmmc[CA].rd_posi[ii]/(float)rdpos_total[ii])*100));}
            else if(kk==CT) {fprintf(fp6,"C>T	%d	%f\n",ii+1,((OA->rmmc[CT].rd_posi[ii]/(float)rdpos_total[ii])*100));}
            else if(kk==CG) {fprintf(fp6,"C>G	%d	%f\n",ii+1,((OA->rmmc[CG].rd_posi[ii]/(float)rdpos_total[ii])*100));}
            else if(kk==GA) {fprintf(fp6,"G>A	%d	%f\n",ii+1,((OA->rmmc[GA].rd_posi[ii]/(float)rdpos_total[ii])*100));}
            else if(kk==GT) {fprintf(fp6,"G>T	%d	%f\n",ii+1,((OA->rmmc[GT].rd_posi[ii]/(float)rdpos_total[ii])*100));}
            else if(kk==GC) {fprintf(fp6,"G>C	%d	%f\n",ii+1,((OA->rmmc[GC].rd_posi[ii]/(float)rdpos_total[ii])*100));}
         }
        
      }

      fprintf(fp7,"Pos	pct_of_reads\n");
      for(ii=0;ii<rLen;ii++) {fprintf(fp7,"%d	%f\n",ii+1,((rdpos_total[ii]/(float)OA->no_indclp)*100));}
      fprintf(fp8,"Pos	Avg_qv\n");
      for(ii=0;ii<rLen;ii++) {fprintf(fp8,"%d	%f\n",ii+1,(OA->rmmQV.rd_posi[ii]/(float)rdpos_total[ii]));}
      
      fprintf(fp9,"Change	Flanks	Pct\n");
      for(kk=0;kk<12;kk++)
      {
         
         if(kk==AT) 
         {
            fprintf(fp9,"A>T	A_A	%f\n",(float)(OA->MMC_flnks[AT][12])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	A_T	%f\n",(float)(OA->MMC_flnks[AT][0])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	A_C	%f\n",(float)(OA->MMC_flnks[AT][1])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	A_G	%f\n",(float)(OA->MMC_flnks[AT][2])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	T_A	%f\n",(float)(OA->MMC_flnks[AT][3])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	T_T	%f\n",(float)(OA->MMC_flnks[AT][13])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	T_C	%f\n",(float)(OA->MMC_flnks[AT][4])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	T_G	%f\n",(float)(OA->MMC_flnks[AT][5])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	C_A	%f\n",(float)(OA->MMC_flnks[AT][6])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	C_T	%f\n",(float)(OA->MMC_flnks[AT][7])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	C_C	%f\n",(float)(OA->MMC_flnks[AT][15])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	C_G	%f\n",(float)(OA->MMC_flnks[AT][8])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	G_A	%f\n",(float)(OA->MMC_flnks[AT][9])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	G_T	%f\n",(float)(OA->MMC_flnks[AT][10])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	G_C	%f\n",(float)(OA->MMC_flnks[AT][11])/(float)OA->MMC[AT]);
            fprintf(fp9,"A>T	G_G	%f\n",(float)(OA->MMC_flnks[AT][14])/(float)OA->MMC[AT]); }
         else if(kk==AC) {
            fprintf(fp9,"A>C	A_A	%f\n",(float)(OA->MMC_flnks[AC][12])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	A_T	%f\n",(float)(OA->MMC_flnks[AC][0])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	A_C	%f\n",(float)(OA->MMC_flnks[AC][1])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	A_G	%f\n",(float)(OA->MMC_flnks[AC][2])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	T_A	%f\n",(float)(OA->MMC_flnks[AC][3])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	T_T	%f\n",(float)(OA->MMC_flnks[AC][13])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	T_C	%f\n",(float)(OA->MMC_flnks[AC][4])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	T_G	%f\n",(float)(OA->MMC_flnks[AC][5])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	C_A	%f\n",(float)(OA->MMC_flnks[AC][6])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	C_T	%f\n",(float)(OA->MMC_flnks[AC][7])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	C_C	%f\n",(float)(OA->MMC_flnks[AC][15])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	C_G	%f\n",(float)(OA->MMC_flnks[AC][8])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	G_A	%f\n",(float)(OA->MMC_flnks[AC][9])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	G_T	%f\n",(float)(OA->MMC_flnks[AC][10])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	G_C	%f\n",(float)(OA->MMC_flnks[AC][11])/(float)OA->MMC[AC]);
            fprintf(fp9,"A>C	G_G	%f\n",(float)(OA->MMC_flnks[AC][14])/(float)OA->MMC[AC]); }
         else if(kk==AG) {
            fprintf(fp9,"A>G	A_A	%f\n",(float)(OA->MMC_flnks[AG][12])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	A_T	%f\n",(float)(OA->MMC_flnks[AG][0])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	A_C	%f\n",(float)(OA->MMC_flnks[AG][1])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	A_G	%f\n",(float)(OA->MMC_flnks[AG][2])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	T_A	%f\n",(float)(OA->MMC_flnks[AG][3])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	T_T	%f\n",(float)(OA->MMC_flnks[AG][13])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	T_C	%f\n",(float)(OA->MMC_flnks[AG][4])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	T_G	%f\n",(float)(OA->MMC_flnks[AG][5])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	C_A	%f\n",(float)(OA->MMC_flnks[AG][6])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	C_T	%f\n",(float)(OA->MMC_flnks[AG][7])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	C_C	%f\n",(float)(OA->MMC_flnks[AG][15])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	C_G	%f\n",(float)(OA->MMC_flnks[AG][8])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	G_A	%f\n",(float)(OA->MMC_flnks[AG][9])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	G_T	%f\n",(float)(OA->MMC_flnks[AG][10])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	G_C	%f\n",(float)(OA->MMC_flnks[AG][11])/(float)OA->MMC[AG]);
            fprintf(fp9,"A>G	G_G	%f\n",(float)(OA->MMC_flnks[AG][14])/(float)OA->MMC[AG]); }
         else if(kk==TA) {
            fprintf(fp9,"T>A	A_A	%f\n",(float)(OA->MMC_flnks[TA][12])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	A_T	%f\n",(float)(OA->MMC_flnks[TA][0])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	A_C	%f\n",(float)(OA->MMC_flnks[TA][1])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	A_G	%f\n",(float)(OA->MMC_flnks[TA][2])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	T_A	%f\n",(float)(OA->MMC_flnks[TA][3])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	T_T	%f\n",(float)(OA->MMC_flnks[TA][13])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	T_C	%f\n",(float)(OA->MMC_flnks[TA][4])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	T_G	%f\n",(float)(OA->MMC_flnks[TA][5])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	C_A	%f\n",(float)(OA->MMC_flnks[TA][6])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	C_T	%f\n",(float)(OA->MMC_flnks[TA][7])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	C_C	%f\n",(float)(OA->MMC_flnks[TA][15])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	C_G	%f\n",(float)(OA->MMC_flnks[TA][8])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	G_A	%f\n",(float)(OA->MMC_flnks[TA][9])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	G_T	%f\n",(float)(OA->MMC_flnks[TA][10])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	G_C	%f\n",(float)(OA->MMC_flnks[TA][11])/(float)OA->MMC[TA]);
            fprintf(fp9,"T>A	G_G	%f\n",(float)(OA->MMC_flnks[TA][14])/(float)OA->MMC[TA]); }
         else if(kk==TC) {
            fprintf(fp9,"T>C	A_A	%f\n",(float)(OA->MMC_flnks[TC][12])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	A_T	%f\n",(float)(OA->MMC_flnks[TC][0])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	A_C	%f\n",(float)(OA->MMC_flnks[TC][1])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	A_G	%f\n",(float)(OA->MMC_flnks[TC][2])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	T_A	%f\n",(float)(OA->MMC_flnks[TC][3])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	T_T	%f\n",(float)(OA->MMC_flnks[TC][13])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	T_C	%f\n",(float)(OA->MMC_flnks[TC][4])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	T_G	%f\n",(float)(OA->MMC_flnks[TC][5])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	C_A	%f\n",(float)(OA->MMC_flnks[TC][6])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	C_T	%f\n",(float)(OA->MMC_flnks[TC][7])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	C_C	%f\n",(float)(OA->MMC_flnks[TC][15])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	C_G	%f\n",(float)(OA->MMC_flnks[TC][8])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	G_A	%f\n",(float)(OA->MMC_flnks[TC][9])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	G_T	%f\n",(float)(OA->MMC_flnks[TC][10])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	G_C	%f\n",(float)(OA->MMC_flnks[TC][11])/(float)OA->MMC[TC]);
            fprintf(fp9,"T>C	G_G	%f\n",(float)(OA->MMC_flnks[TC][14])/(float)OA->MMC[TC]); }
         else if(kk==TG) {
            fprintf(fp9,"T>G	A_A	%f\n",(float)(OA->MMC_flnks[TG][12])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	A_T	%f\n",(float)(OA->MMC_flnks[TG][0])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	A_C	%f\n",(float)(OA->MMC_flnks[TG][1])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	A_G	%f\n",(float)(OA->MMC_flnks[TG][2])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	T_A	%f\n",(float)(OA->MMC_flnks[TG][3])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	T_T	%f\n",(float)(OA->MMC_flnks[TG][13])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	T_C	%f\n",(float)(OA->MMC_flnks[TG][4])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	T_G	%f\n",(float)(OA->MMC_flnks[TG][5])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	C_A	%f\n",(float)(OA->MMC_flnks[TG][6])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	C_T	%f\n",(float)(OA->MMC_flnks[TG][7])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	C_C	%f\n",(float)(OA->MMC_flnks[TG][15])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	C_G	%f\n",(float)(OA->MMC_flnks[TG][8])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	G_A	%f\n",(float)(OA->MMC_flnks[TG][9])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	G_T	%f\n",(float)(OA->MMC_flnks[TG][10])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	G_C	%f\n",(float)(OA->MMC_flnks[TG][11])/(float)OA->MMC[TG]);
            fprintf(fp9,"T>G	G_G	%f\n",(float)(OA->MMC_flnks[TG][14])/(float)OA->MMC[TG]); }
         else if(kk==CA) {
            fprintf(fp9,"C>A	A_A	%f\n",(float)(OA->MMC_flnks[CA][12])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	A_T	%f\n",(float)(OA->MMC_flnks[CA][0])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	A_C	%f\n",(float)(OA->MMC_flnks[CA][1])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	A_G	%f\n",(float)(OA->MMC_flnks[CA][2])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	T_A	%f\n",(float)(OA->MMC_flnks[CA][3])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	T_T	%f\n",(float)(OA->MMC_flnks[CA][13])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	T_C	%f\n",(float)(OA->MMC_flnks[CA][4])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	T_G	%f\n",(float)(OA->MMC_flnks[CA][5])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	C_A	%f\n",(float)(OA->MMC_flnks[CA][6])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	C_T	%f\n",(float)(OA->MMC_flnks[CA][7])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	C_C	%f\n",(float)(OA->MMC_flnks[CA][15])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	C_G	%f\n",(float)(OA->MMC_flnks[CA][8])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	G_A	%f\n",(float)(OA->MMC_flnks[CA][9])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	G_T	%f\n",(float)(OA->MMC_flnks[CA][10])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	G_C	%f\n",(float)(OA->MMC_flnks[CA][11])/(float)OA->MMC[CA]);
            fprintf(fp9,"C>A	G_G	%f\n",(float)(OA->MMC_flnks[CA][14])/(float)OA->MMC[CA]); }
         else if(kk==CT) {
            fprintf(fp9,"C>T	A_A	%f\n",(float)(OA->MMC_flnks[CT][12])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	A_T	%f\n",(float)(OA->MMC_flnks[CT][0])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	A_C	%f\n",(float)(OA->MMC_flnks[CT][1])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	A_G	%f\n",(float)(OA->MMC_flnks[CT][2])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	T_A	%f\n",(float)(OA->MMC_flnks[CT][3])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	T_T	%f\n",(float)(OA->MMC_flnks[CT][13])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	T_C	%f\n",(float)(OA->MMC_flnks[CT][4])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	T_G	%f\n",(float)(OA->MMC_flnks[CT][5])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	C_A	%f\n",(float)(OA->MMC_flnks[CT][6])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	C_T	%f\n",(float)(OA->MMC_flnks[CT][7])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	C_C	%f\n",(float)(OA->MMC_flnks[CT][15])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	C_G	%f\n",(float)(OA->MMC_flnks[CT][8])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	G_A	%f\n",(float)(OA->MMC_flnks[CT][9])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	G_T	%f\n",(float)(OA->MMC_flnks[CT][10])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	G_C	%f\n",(float)(OA->MMC_flnks[CT][11])/(float)OA->MMC[CT]);
            fprintf(fp9,"C>T	G_G	%f\n",(float)(OA->MMC_flnks[CT][14])/(float)OA->MMC[CT]); }
         else if(kk==CG) {
            fprintf(fp9,"C>G	A_A	%f\n",(float)(OA->MMC_flnks[CG][12])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	A_T	%f\n",(float)(OA->MMC_flnks[CG][0])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	A_C	%f\n",(float)(OA->MMC_flnks[CG][1])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	A_G	%f\n",(float)(OA->MMC_flnks[CG][2])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	T_A	%f\n",(float)(OA->MMC_flnks[CG][3])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	T_T	%f\n",(float)(OA->MMC_flnks[CG][13])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	T_C	%f\n",(float)(OA->MMC_flnks[CG][4])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	T_G	%f\n",(float)(OA->MMC_flnks[CG][5])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	C_A	%f\n",(float)(OA->MMC_flnks[CG][6])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	C_T	%f\n",(float)(OA->MMC_flnks[CG][7])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	C_C	%f\n",(float)(OA->MMC_flnks[CG][15])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	C_G	%f\n",(float)(OA->MMC_flnks[CG][8])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	G_A	%f\n",(float)(OA->MMC_flnks[CG][9])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	G_T	%f\n",(float)(OA->MMC_flnks[CG][10])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	G_C	%f\n",(float)(OA->MMC_flnks[CG][11])/(float)OA->MMC[CG]);
            fprintf(fp9,"C>G	G_G	%f\n",(float)(OA->MMC_flnks[CG][14])/(float)OA->MMC[CG]); }
         else if(kk==GA) {
            fprintf(fp9,"G>A	A_A	%f\n",(float)(OA->MMC_flnks[GA][12])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	A_T	%f\n",(float)(OA->MMC_flnks[GA][0])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	A_C	%f\n",(float)(OA->MMC_flnks[GA][1])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	A_G	%f\n",(float)(OA->MMC_flnks[GA][2])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	T_A	%f\n",(float)(OA->MMC_flnks[GA][3])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	T_T	%f\n",(float)(OA->MMC_flnks[GA][13])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	T_C	%f\n",(float)(OA->MMC_flnks[GA][4])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	T_G	%f\n",(float)(OA->MMC_flnks[GA][5])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	C_A	%f\n",(float)(OA->MMC_flnks[GA][6])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	C_T	%f\n",(float)(OA->MMC_flnks[GA][7])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	C_C	%f\n",(float)(OA->MMC_flnks[GA][15])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	C_G	%f\n",(float)(OA->MMC_flnks[GA][8])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	G_A	%f\n",(float)(OA->MMC_flnks[GA][9])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	G_T	%f\n",(float)(OA->MMC_flnks[GA][10])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	G_C	%f\n",(float)(OA->MMC_flnks[GA][11])/(float)OA->MMC[GA]);
            fprintf(fp9,"G>A	G_G	%f\n",(float)(OA->MMC_flnks[GA][14])/(float)OA->MMC[GA]); }
         else if(kk==GT) {
            fprintf(fp9,"G>T	A_A	%f\n",(float)(OA->MMC_flnks[GT][12])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	A_T	%f\n",(float)(OA->MMC_flnks[GT][0])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	A_C	%f\n",(float)(OA->MMC_flnks[GT][1])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	A_G	%f\n",(float)(OA->MMC_flnks[GT][2])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	T_A	%f\n",(float)(OA->MMC_flnks[GT][3])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	T_T	%f\n",(float)(OA->MMC_flnks[GT][13])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	T_C	%f\n",(float)(OA->MMC_flnks[GT][4])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	T_G	%f\n",(float)(OA->MMC_flnks[GT][5])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	C_A	%f\n",(float)(OA->MMC_flnks[GT][6])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	C_T	%f\n",(float)(OA->MMC_flnks[GT][7])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	C_C	%f\n",(float)(OA->MMC_flnks[GT][15])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	C_G	%f\n",(float)(OA->MMC_flnks[GT][8])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	G_A	%f\n",(float)(OA->MMC_flnks[GT][9])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	G_T	%f\n",(float)(OA->MMC_flnks[GT][10])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	G_C	%f\n",(float)(OA->MMC_flnks[GT][11])/(float)OA->MMC[GT]);
            fprintf(fp9,"G>T	G_G	%f\n",(float)(OA->MMC_flnks[GT][14])/(float)OA->MMC[GT]); }
         else if(kk==GC) {
            fprintf(fp9,"G>C	A_A	%f\n",(float)(OA->MMC_flnks[GC][12])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	A_T	%f\n",(float)(OA->MMC_flnks[GC][0])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	A_C	%f\n",(float)(OA->MMC_flnks[GC][1])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	A_G	%f\n",(float)(OA->MMC_flnks[GC][2])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	T_A	%f\n",(float)(OA->MMC_flnks[GC][3])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	T_T	%f\n",(float)(OA->MMC_flnks[GC][13])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	T_C	%f\n",(float)(OA->MMC_flnks[GC][4])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	T_G	%f\n",(float)(OA->MMC_flnks[GC][5])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	C_A	%f\n",(float)(OA->MMC_flnks[GC][6])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	C_T	%f\n",(float)(OA->MMC_flnks[GC][7])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	C_C	%f\n",(float)(OA->MMC_flnks[GC][15])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	C_G	%f\n",(float)(OA->MMC_flnks[GC][8])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	G_A	%f\n",(float)(OA->MMC_flnks[GC][9])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	G_T	%f\n",(float)(OA->MMC_flnks[GC][10])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	G_C	%f\n",(float)(OA->MMC_flnks[GC][11])/(float)OA->MMC[GC]);
            fprintf(fp9,"G>C	G_G	%f\n",(float)(OA->MMC_flnks[GC][14])/(float)OA->MMC[GC]); }
      }
      
      fprintf(fp10,"Length	Flanks	Pct\n");
      for(kk=0;kk<5;kk++)
      {
         
         if(kk==4) 
         {
            fprintf(fp10,">=-5	A_A	%f\n",(float)(OA->tDel_flnks[kk][12])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	A_T	%f\n",(float)(OA->tDel_flnks[kk][0])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	A_C	%f\n",(float)(OA->tDel_flnks[kk][1])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	A_G	%f\n",(float)(OA->tDel_flnks[kk][2])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	T_A	%f\n",(float)(OA->tDel_flnks[kk][3])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	T_T	%f\n",(float)(OA->tDel_flnks[kk][13])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	T_C	%f\n",(float)(OA->tDel_flnks[kk][4])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	T_G	%f\n",(float)(OA->tDel_flnks[kk][5])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	C_A	%f\n",(float)(OA->tDel_flnks[kk][6])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	C_T	%f\n",(float)(OA->tDel_flnks[kk][7])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	C_C	%f\n",(float)(OA->tDel_flnks[kk][15])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	C_G	%f\n",(float)(OA->tDel_flnks[kk][8])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	G_A	%f\n",(float)(OA->tDel_flnks[kk][9])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	G_T	%f\n",(float)(OA->tDel_flnks[kk][10])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	G_C	%f\n",(float)(OA->tDel_flnks[kk][11])/(float)(OA->tDel[kk]));
            fprintf(fp10,">=-5	G_G	%f\n",(float)(OA->tDel_flnks[kk][14])/(float)(OA->tDel[kk])); }
         else 
         {
            fprintf(fp10,"-%d	A_A	%f\n", kk+1, (float)(OA->tDel_flnks[kk][12])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	A_T	%f\n", kk+1, (float)(OA->tDel_flnks[kk][0])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	A_C	%f\n", kk+1, (float)(OA->tDel_flnks[kk][1])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	A_G	%f\n", kk+1, (float)(OA->tDel_flnks[kk][2])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	T_A	%f\n", kk+1, (float)(OA->tDel_flnks[kk][3])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	T_T	%f\n", kk+1, (float)(OA->tDel_flnks[kk][13])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	T_C	%f\n", kk+1, (float)(OA->tDel_flnks[kk][4])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	T_G	%f\n", kk+1, (float)(OA->tDel_flnks[kk][5])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	C_A	%f\n", kk+1, (float)(OA->tDel_flnks[kk][6])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	C_T	%f\n", kk+1, (float)(OA->tDel_flnks[kk][7])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	C_C	%f\n", kk+1, (float)(OA->tDel_flnks[kk][15])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	C_G	%f\n", kk+1, (float)(OA->tDel_flnks[kk][8])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	G_A	%f\n", kk+1, (float)(OA->tDel_flnks[kk][9])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	G_T	%f\n", kk+1, (float)(OA->tDel_flnks[kk][10])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	G_C	%f\n", kk+1, (float)(OA->tDel_flnks[kk][11])/(float)(OA->tDel[kk]));
            fprintf(fp10,"-%d	G_G	%f\n", kk+1, (float)(OA->tDel_flnks[kk][14])/(float)(OA->tDel[kk])); } 
      }
            fprintf(fp10,"0	A_A	0\n");
            fprintf(fp10,"0	A_T	0\n");
            fprintf(fp10,"0	A_C	0\n");
            fprintf(fp10,"0	A_G	0\n");
            fprintf(fp10,"0	T_A	0\n");
            fprintf(fp10,"0	T_T	0\n");
            fprintf(fp10,"0	T_C	0\n");
            fprintf(fp10,"0	T_G	0\n");
            fprintf(fp10,"0	C_A	0\n");
            fprintf(fp10,"0	C_T	0\n");
            fprintf(fp10,"0	C_C	0\n");
            fprintf(fp10,"0	C_G	0\n");
            fprintf(fp10,"0	G_A	0\n");
            fprintf(fp10,"0	G_T	0\n");
            fprintf(fp10,"0	G_C	0\n");
            fprintf(fp10,"0	G_G	0\n");
      
      for(kk=0;kk<5;kk++)
      {
         
         if(kk==4) 
         {
            fprintf(fp10,">=+5	A_A	%f\n",(float)(OA->tIns_flnks[kk][12])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	A_T	%f\n",(float)(OA->tIns_flnks[kk][0])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	A_C	%f\n",(float)(OA->tIns_flnks[kk][1])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	A_G	%f\n",(float)(OA->tIns_flnks[kk][2])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	T_A	%f\n",(float)(OA->tIns_flnks[kk][3])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	T_T	%f\n",(float)(OA->tIns_flnks[kk][13])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	T_C	%f\n",(float)(OA->tIns_flnks[kk][4])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	T_G	%f\n",(float)(OA->tIns_flnks[kk][5])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	C_A	%f\n",(float)(OA->tIns_flnks[kk][6])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	C_T	%f\n",(float)(OA->tIns_flnks[kk][7])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	C_C	%f\n",(float)(OA->tIns_flnks[kk][15])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	C_G	%f\n",(float)(OA->tIns_flnks[kk][8])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	G_A	%f\n",(float)(OA->tIns_flnks[kk][9])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	G_T	%f\n",(float)(OA->tIns_flnks[kk][10])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	G_C	%f\n",(float)(OA->tIns_flnks[kk][11])/(float)(OA->tIns[kk]));
            fprintf(fp10,">=+5	G_G	%f\n",(float)(OA->tIns_flnks[kk][14])/(float)(OA->tIns[kk])); }
         else 
         {
            fprintf(fp10,"+%d	A_A	%f\n", kk+1, (float)(OA->tIns_flnks[kk][12])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	A_T	%f\n", kk+1, (float)(OA->tIns_flnks[kk][0])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	A_C	%f\n", kk+1, (float)(OA->tIns_flnks[kk][1])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	A_G	%f\n", kk+1, (float)(OA->tIns_flnks[kk][2])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	T_A	%f\n", kk+1, (float)(OA->tIns_flnks[kk][3])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	T_T	%f\n", kk+1, (float)(OA->tIns_flnks[kk][13])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	T_C	%f\n", kk+1, (float)(OA->tIns_flnks[kk][4])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	T_G	%f\n", kk+1, (float)(OA->tIns_flnks[kk][5])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	C_A	%f\n", kk+1, (float)(OA->tIns_flnks[kk][6])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	C_T	%f\n", kk+1, (float)(OA->tIns_flnks[kk][7])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	C_C	%f\n", kk+1, (float)(OA->tIns_flnks[kk][15])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	C_G	%f\n", kk+1, (float)(OA->tIns_flnks[kk][8])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	G_A	%f\n", kk+1, (float)(OA->tIns_flnks[kk][9])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	G_T	%f\n", kk+1, (float)(OA->tIns_flnks[kk][10])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	G_C	%f\n", kk+1, (float)(OA->tIns_flnks[kk][11])/(float)(OA->tIns[kk]));
            fprintf(fp10,"+%d	G_G	%f\n", kk+1, (float)(OA->tIns_flnks[kk][14])/(float)(OA->tIns[kk])); } 
      }
      
      
      
      fclose(fp1);
      fclose(fp1c);
      fclose(fp2);
      fclose(fp3);
      fclose(fp4);
      fclose(fp5);
      fclose(fp6);
      fclose(fp7);
      fclose(fp8);
      fclose(fp9);
      fclose(fp10);
      free(outfn);
     
}


 
 void PrintGenrlLogFiles_u(log_v *T, int rLen, long int trgtsize, char *outdir, char *prefix, paramtr *prmtr)
 {
        FILE *gfp1;
        int ii;
        long int atcg_cnt=0;
        char *outfn=0;
        outfn = calloc(strlen(outdir) + 500, 1);
        sprintf(outfn, "%s/%s_mapping_summary.log", outdir, prefix);
        gfp1=fopen(outfn,"w");
        //rd_cnt=T->glog.tot_rd_cnt-(T->glog.scndry_cnt + T->glog.scndry_cnt);
        fprintf(gfp1,"\nmapinsights bamqc report (version-1.0)\n");
        fprintf(gfp1,"======================================\n\n");
        fprintf(gfp1,"command = mapinsights %s\n\n", prmtr->comnd);
        fprintf(gfp1,"Mapping Statistics\n");
        fprintf(gfp1,"-------------------\n");
        if(!strcmp(prefix, "Overall"))
        {
        fprintf(gfp1,"Total number of reads = %ld\n",T->glog.all_rd_cnt);
        fprintf(gfp1,"Total number of mapped reads = %ld (%.2f%%)\n",(T->glog.all_rd_cnt-T->glog.unmap), ((float)(T->glog.all_rd_cnt-T->glog.unmap)/(float)T->glog.all_rd_cnt)*100);
        fprintf(gfp1,"Total number of unmapped reads = %ld (%.2f%%)\n",T->glog.unmap, ((float)(T->glog.unmap)/(float)T->glog.all_rd_cnt)*100);
        }
        fprintf(gfp1,"\n");
        if(trgtsize) {fprintf(gfp1,"Total number of mapped reads in targeted regions = %ld \n",T->glog.tot_rd_cnt);}
        //else {fprintf(gfp1,"Total number of mapped reads= %ld (%.2f)\n",T->glog.tot_rd_cnt, ((float)(T->glog.tot_rd_cnt)/(float)T->glog.tot_rd_cnt)*100);}
        fprintf(gfp1,"Mapped-read1 = %ld (%.2f%%)\n", T->glog.r1_cnt, ((float)(T->glog.r1_cnt)/(float)T->glog.tot_rd_cnt)*100);
	fprintf(gfp1,"Mapped-read2 = %ld (%.2f%%)\n",T->glog.r2_cnt, ((float)(T->glog.r2_cnt)/(float)T->glog.tot_rd_cnt)*100);
	fprintf(gfp1,"Mapped-forward = %ld (%.2f%%)\n",(T->glog.tot_rd_cnt-T->glog.rev_cnt), ((float)((T->glog.tot_rd_cnt-T->glog.rev_cnt))/(float)T->glog.tot_rd_cnt)*100);
	fprintf(gfp1,"Mapped-reverse = %ld (%.2f%%)\n",T->glog.rev_cnt, ((float)(T->glog.rev_cnt)/(float)T->glog.tot_rd_cnt)*100);
	fprintf(gfp1,"QC-failed = %ld (%.2f%%)\n",T->glog.qc_fail, ((float)(T->glog.qc_fail)/(float)T->glog.tot_rd_cnt)*100);
	fprintf(gfp1,"Mapped-pair = %ld (%.2f%%)\n",T->glog.map_pr, ((float)(T->glog.map_pr)/(float)T->glog.tot_rd_cnt)*100);	
	fprintf(gfp1,"Mapped-properpair = %ld (%.2f%%)\n",T->glog.propr_pr, ((float)(T->glog.propr_pr)/(float)T->glog.tot_rd_cnt)*100);
	fprintf(gfp1,"Secondary-alignments = %ld (%.2f%%)\n",T->glog.scndry_cnt, ((float)(T->glog.scndry_cnt)/(float)T->glog.tot_rd_cnt)*100);
	fprintf(gfp1,"Supplementary-alignments = %ld (%.2f%%)\n",T->glog.supply_cnt, ((float)(T->glog.supply_cnt)/(float)T->glog.tot_rd_cnt)*100);
	fprintf(gfp1,"Strand ratio(F:R) = %.2f:%.2f\n",(float)(T->glog.tot_rd_cnt-T->glog.rev_cnt)/(float)T->glog.tot_rd_cnt,(float)T->glog.rev_cnt/(float)T->glog.tot_rd_cnt);
	
	fprintf(gfp1,"\n");
	fprintf(gfp1,"Softclip-events = %ld (%.2f%%)\n", T->midc.sftclp, ((float)(T->midc.sftclp)/(float)(T->midc.sftclp+T->midc.hrdclp))*100);
	fprintf(gfp1,"Softclipped basecounts = %ld\n", T->midc.sftclp_bs);
	fprintf(gfp1,"Hardclip-events = %ld (%.2f%%)\n", T->midc.hrdclp, ((float)(T->midc.hrdclp)/(float)(T->midc.sftclp+T->midc.hrdclp))*100);
	fprintf(gfp1,"Hardclipped basecounts = %ld\n", T->midc.hrdclp_bs);
	
	fprintf(gfp1,"\n");
	fprintf(gfp1,"Nextera-Transposase-Sequence = %ld (%.2f%%)\n",T->adptr.NTS, ((float)(T->adptr.NTS)/(float)T->glog.tot_rd_cnt)*100);
	fprintf(gfp1,"Illumina-PCR-primer = %ld (%.2f%%)\n",T->adptr.IPCRprm, ((float)(T->adptr.IPCRprm)/(float)T->glog.tot_rd_cnt)*100);
	fprintf(gfp1,"Illumina-Adapter = %ld (%.2f%%)\n",T->adptr.IADptr, ((float)(T->adptr.IADptr)/(float)T->glog.tot_rd_cnt)*100);
	
	fprintf(gfp1,"\n");
	fprintf(gfp1,"Aligned-nucleotides = %ld\n", (T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt));
	fprintf(gfp1,"A = %ld (%.2f%%)\n",T->glog.AAcnt ,((double)(T->glog.AAcnt)/(double)(T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt)*100));
	fprintf(gfp1,"T = %ld (%.2f%%)\n",T->glog.TTcnt ,((double)(T->glog.TTcnt)/(double)(T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt)*100));
	fprintf(gfp1,"G = %ld (%.2f%%)\n",T->glog.GGcnt ,((double)(T->glog.GGcnt)/(double)(T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt)*100));
	fprintf(gfp1,"C = %ld (%.2f%%)\n",T->glog.CCcnt ,((double)(T->glog.CCcnt)/(double)(T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt)*100));
	fprintf(gfp1,"N = %ld (%.2f%%)\n",T->glog.NNcnt ,((double)(T->glog.NNcnt)/(double)(T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt)*100));
	fprintf(gfp1,"\n");
	fprintf(gfp1,"GC %% = %f\n",((double)(T->glog.GGcnt+T->glog.CCcnt)/(double)(T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt)*100));
        fprintf(gfp1,"\n");
        
        
	fprintf(gfp1,"Single nucleotide matchs = %ld\n", (long int)(T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt)-(T->midc.basemmc));
	fprintf(gfp1,"Single nucleotide mismatchs = %ld\n", T->midc.basemmc);
	fprintf(gfp1,"Insertion-events = %ld\n", T->midc.ins_cnt);
	fprintf(gfp1,"Inserted basecount = %ld\n", T->midc.ins_cnt_bs);
	fprintf(gfp1,"Homopolymer Insertation = %ld (%.2f%%)\n", T->midc.polyIns, ((double)(T->midc.polyIns)/(double)(T->midc.ins_cnt)*100));
	fprintf(gfp1,"Deletion-events = %ld\n", T->midc.del_cnt);
	fprintf(gfp1,"Deleted basecounts = %ld\n", T->midc.del_cnt_bs);
	fprintf(gfp1,"Homopolymer Deletion = %ld (%.2f%%)\n", T->midc.polyDel, ((double)(T->midc.polyDel)/(double)(T->midc.del_cnt)*100));
	fprintf(gfp1,"Mismatch rate = %f\n", (float)(T->midc.basemmc)/(float)(T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt));
	fprintf(gfp1,"Insertion rate = %f\n", (float)(T->midc.ins_cnt)/(float)(T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt));
	fprintf(gfp1,"Deletion rate = %f\n", (float)(T->midc.del_cnt)/(float)(T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt));
	if(!strcmp(prefix, "Overall"))
        {
	fprintf(gfp1,"\n");
	fprintf(gfp1,"Mapped-diffChr = %ld (%.2f%%)\n",T->glog.map_difchr, ((float)(T->glog.map_difchr)/(float)T->glog.tot_rd_cnt)*100);
	fprintf(gfp1,"Insertsize>=1K = %ld (%.2f%%)\n",T->glog.Isize_ge1K, ((float)(T->glog.Isize_ge1K)/(float)T->glog.tot_rd_cnt)*100);
	fprintf(gfp1,"\n");
	fprintf(gfp1,"Mapped-read orientation RR:FF:Innerward:Outward = %ld (%.2f%%):%ld (%.2f%%):%ld (%.2f%%):%ld (%.2f%%)\n",T->glog.RR, ((float)(T->glog.RR)/(float)T->glog.tot_rd_cnt)*100,T->glog.FF, ((float)(T->glog.FF)/(float)T->glog.tot_rd_cnt)*100,(T->glog.FR-T->glog.OO), ((float)(T->glog.FR-T->glog.OO)/(float)T->glog.tot_rd_cnt)*100,T->glog.OO, ((float)(T->glog.OO)/(float)T->glog.tot_rd_cnt)*100);
	
	fprintf(gfp1,"\n");
	fprintf(gfp1,"Mean depth of coverages = %f\n",(float)(T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt)/(float)trgtsize);
	fprintf(gfp1,"DuplicateMarked rate = %ld (%.2f%%)\n", T->glog.DupCnt, ((float)(T->glog.DupCnt)/(float)T->glog.tot_rd_cnt)*100);
        }
	fclose(gfp1);
	
      if(!strcmp(prefix, "Overall"))
      {
      *outfn=0;
      sprintf(outfn, "%s/Overall_Basequality_dist.txt", outdir);
      gfp1=fopen(outfn,"w");
      fprintf(gfp1,"BaseQ	Basecounts\n");
      for(ii=0;ii<=MAX_BASEQ;ii++) {fprintf(gfp1,"%d	%ld\n",ii,T->glog.BaseQ_dist[ii]);}
      fclose(gfp1);
      //printf("MapQ	No_of_reads\n");
      //printf("0	%d	%f\n", T->glog.MapQ_bin_dist[1], ((float)T->glog.MapQ_bin_dist[1]/(float)T->glog.tot_rd_cnt)*100);
      //printf("1-10	%d	%f\n", T->glog.MapQ_bin_dist[2], ((float)T->glog.MapQ_bin_dist[2]/(float)T->glog.tot_rd_cnt)*100);
      //printf("11-20	%d	%f\n", T->glog.MapQ_bin_dist[3], ((float)T->glog.MapQ_bin_dist[3]/(float)T->glog.tot_rd_cnt)*100);
      //printf("21-30	%d	%f\n", T->glog.MapQ_bin_dist[4], ((float)T->glog.MapQ_bin_dist[4]/(float)T->glog.tot_rd_cnt)*100);
      //printf("31-40	%d	%f\n", T->glog.MapQ_bin_dist[5], ((float)T->glog.MapQ_bin_dist[5]/(float)T->glog.tot_rd_cnt)*100);
      //printf("41-50	%d	%f\n", T->glog.MapQ_bin_dist[6], ((float)T->glog.MapQ_bin_dist[6]/(float)T->glog.tot_rd_cnt)*100);
      //printf(">50	%d	%f\n", T->glog.MapQ_bin_dist[7], ((float)T->glog.MapQ_bin_dist[7]/(float)T->glog.tot_rd_cnt)*100);
      
      *outfn=0;
      sprintf(outfn, "%s/Overall_Mappingquality_bin_dist.txt", outdir);
      gfp1=fopen(outfn,"w");
      fprintf(gfp1,"MapQ	Pct_reads\n");
      fprintf(gfp1,"0	%f\n", ((float)T->glog.MapQ_bin_dist[1]/(float)T->glog.tot_rd_cnt)*100);
      fprintf(gfp1,"1-10	%f\n", ((float)T->glog.MapQ_bin_dist[2]/(float)T->glog.tot_rd_cnt)*100);
      fprintf(gfp1,"11-20	%f\n", ((float)T->glog.MapQ_bin_dist[3]/(float)T->glog.tot_rd_cnt)*100);
      fprintf(gfp1,"21-30	%f\n", ((float)T->glog.MapQ_bin_dist[4]/(float)T->glog.tot_rd_cnt)*100);
      fprintf(gfp1,"31-40	%f\n", ((float)T->glog.MapQ_bin_dist[5]/(float)T->glog.tot_rd_cnt)*100);
      fprintf(gfp1,"41-50	%f\n", ((float)T->glog.MapQ_bin_dist[6]/(float)T->glog.tot_rd_cnt)*100);
      fprintf(gfp1,"51-more	%f\n", ((float)T->glog.MapQ_bin_dist[7]/(float)T->glog.tot_rd_cnt)*100); 
      fclose(gfp1);

      *outfn=0;
      sprintf(outfn, "%s/Overall_Mappingquality_dist.txt", outdir);
      gfp1=fopen(outfn,"w");
      fprintf(gfp1,"MapQ	No_of_reads\n");
      for(ii=1;ii<=MAX_MAPQ;ii++) {fprintf(gfp1,"%d	%ld\n",ii,T->glog.MapQ_dist[ii]);}
      fclose(gfp1);

      *outfn=0;
      sprintf(outfn, "%s/Overall_Insertsize_dist.txt", outdir);
      gfp1=fopen(outfn,"w");
      fprintf(gfp1,"Isize	No_of_reads\n");
      for(ii=1;ii<MAX_ISIZE;ii++) {fprintf(gfp1,"%d	%ld\n",ii,T->glog.Isize_dist[ii]);}
      fclose(gfp1);
      
      *outfn=0;
      sprintf(outfn, "%s/Overall_GC_dist.txt", outdir);
      gfp1=fopen(outfn,"w");
      fprintf(gfp1,"GCpct	pct_of_reads\n");
      for(ii=0;ii<100;ii=ii+2) {fprintf(gfp1,"%d	%f\n",ii,((float)(T->glog.GC_pct[ii]+T->glog.GC_pct[ii+1])/(float)T->glog.tot_rd_cnt)*100);}

      fprintf(gfp1,"100	%f\n",((float)T->glog.GC_pct[100]/(float)T->glog.tot_rd_cnt)*100);
      fclose(gfp1);
      
      *outfn=0;
      sprintf(outfn, "%s/Overall_Perbase_seq-content.txt", outdir);
      gfp1=fopen(outfn,"w");
      fprintf(gfp1,"Readpos	Base	Pct\n");
      for(ii=0;ii<rLen;ii++) {
      atcg_cnt = T->glog.A_con[ii] + T->glog.T_con[ii] + T->glog.G_con[ii] + T->glog.C_con[ii] + T->glog.N_con[ii];
      if(atcg_cnt == 0) continue;
      fprintf(gfp1,"%d	A	%f\n%d	T	%f\n%d	G	%f\n%d	C	%f\n%d	N	%f\n",ii+1,((float)T->glog.A_con[ii]/(float)atcg_cnt)*100,ii+1,((float)T->glog.T_con[ii]/(float)atcg_cnt)*100,ii+1,((float)T->glog.G_con[ii]/(float)atcg_cnt)*100,ii+1,((float)T->glog.C_con[ii]/(float)atcg_cnt)*100,ii+1,((float)T->glog.N_con[ii]/(float)atcg_cnt)*100);
       }
      //for(ii=0;ii<rLen;ii++) {fprintf(gfp1,"%d	A	%f\n%d	T	%f\n%d	G	%f\n%d	C	%f\n%d	N	%f\n",ii+1,((float)T->glog.A_con[ii]/(float)rd_cnt)*100,ii+1,((float)T->glog.T_con[ii]/(float)rd_cnt)*100,ii+1,((float)T->glog.G_con[ii]/(float)rd_cnt)*100,ii+1,((float)T->glog.C_con[ii]/(float)rd_cnt)*100,ii+1,((float)T->glog.N_con[ii]/(float)rd_cnt)*100);}
      fclose(gfp1);
      
      *outfn=0;
      sprintf(outfn, "%s/Overall_Percycle_avg_baseQ.txt", outdir);
      gfp1=fopen(outfn,"w");
      fprintf(gfp1,"Readpos	Avg_baseQ\n");      
      for(ii=0;ii<rLen;ii++) {fprintf(gfp1,"%d	%f\n",ii+1,(float)T->glog.Percycle_mn_baseQ[ii]/(float)T->glog.tot_rd_cnt);}
      fclose(gfp1);
      
      *outfn=0;
      sprintf(outfn, "%s/Overall_Readlength_dist.txt", outdir);
      gfp1=fopen(outfn,"w");
      fprintf(gfp1,"Length	No_of_reads\n");      
      for(ii=1;ii<=rLen;ii++) {fprintf(gfp1,"%d	%ld\n",ii,T->glog.Rdlen[ii]);}
      fclose(gfp1);

       *outfn=0;
      sprintf(outfn, "%s/Add_log.r", outdir);
      gfp1=fopen(outfn,"w");
      fprintf(gfp1,"isize <- read.table(\"%s/Overall_Insertsize_dist.txt\",header = TRUE, sep = \"\\t\")\n", outdir);

      fprintf(gfp1,"calmean <- function(isize)\n");
      fprintf(gfp1,"with(isize, as.double(sum(as.double(Isize*(No_of_reads/10000))))/as.double(as.double(sum(No_of_reads/10000))))\n");
      fprintf(gfp1,"mean <- paste(\"Mean insertsize =\", round(calmean(isize), 2))\n");
      fprintf(gfp1,"write(mean, file=\"%s/Overall_mapping_summary.log\", append=TRUE)\n", outdir);
   
      fprintf(gfp1,"baseq <- read.table(\"%s/Overall_Basequality_dist.txt\",header = TRUE, sep = \"\\t\")\n", outdir);
      fprintf(gfp1,"calmean <- function(baseq)\n");
      fprintf(gfp1,"with(baseq, as.double(sum(as.double(BaseQ*(Basecounts/100000))))/as.double(as.double(sum(Basecounts/100000))))\n");
      //fprintf(gfp1,"with(baseq, as.double(sum(as.double(BaseQ*Basecounts))/as.numeric(sum(Basecounts)))\n");
      //fprintf(gfp1,"round(calmean(baseq),2)\n");
      fprintf(gfp1,"mean <- paste(\"Mean basequality =\", round(calmean(baseq),2))\n");
      fprintf(gfp1,"write(mean, file=\"%s/Overall_mapping_summary.log\", append=TRUE)\n", outdir);

      fprintf(gfp1,"mapq <- read.table(\"%s/Overall_Mappingquality_dist.txt\",header = TRUE, sep = \"\\t\")\n", outdir);
      fprintf(gfp1,"calmean <- function(mapq)\n");
      fprintf(gfp1,"with(mapq, as.double(sum(as.double(MapQ*(No_of_reads/10000))))/as.double(as.double(sum(No_of_reads/10000))))\n");
      //fprintf(gfp1,"with(mapq, as.numeric(sum(MapQ*No_of_reads))/as.numeric(sum(No_of_reads)))\n");
      fprintf(gfp1,"mean <- paste(\"Mean mapping quality =\" ,round(calmean(mapq),2))\n");
      fprintf(gfp1,"write(mean, file=\"%s/Overall_mapping_summary.log\", append=TRUE)\n", outdir);

      fprintf(gfp1,"rlen <- read.table(\"%s/Overall_Readlength_dist.txt\",header = TRUE, sep = \"\\t\")\n", outdir);
      fprintf(gfp1,"calmean <- function(rlen)\n");
      fprintf(gfp1,"with(rlen, as.double(sum(as.double(Length*(No_of_reads/10000))))/as.double(as.double(sum(No_of_reads/10000))))\n");
      //fprintf(gfp1,"with(rlen, as.numeric(sum(Length*No_of_reads))/as.numeric(sum(No_of_reads)))\n");
      fprintf(gfp1,"mean <- paste(\"Mean read length =\" ,round(calmean(rlen),2))\n");
      fprintf(gfp1,"write(mean, file=\"%s/Overall_mapping_summary.log\", append=TRUE)\n", outdir);
      
      fclose(gfp1);
   
      *outfn=0;
      sprintf(outfn, "R CMD BATCH %s/Add_log.r", outdir);
      system(outfn);
      
      *outfn=0;
      sprintf(outfn, "rm -f %s/Add_log.*", outdir);
      system(outfn);
      system("rm -f Add_log.*");
      //free(outfn);
   }
   free(outfn);
 }
 void draw_plots(char *outdir)
 {
   FILE *gfp1;
   char *outfn=0;
   outfn = calloc(strlen(outdir) + 500, 1);
   sprintf(outfn, "%s/statistical_test.txt", outdir);
   gfp1=fopen(outfn,"w");
   fclose(gfp1);
   *outfn=0;
   outfn = calloc(strlen(outdir) + 500, 1);
   sprintf(outfn, "%s/plotsmake.r", outdir);
   //printf("Outfile : %s\n",outfn);
   gfp1=fopen(outfn,"w");
   fprintf(gfp1,"library (\"ggplot2\")\n");
   
   ////*updated on 19Apr2021 start*///
   fprintf(gfp1,"r1_fs <- read.table(\"%s/Read1_fs_SN_MM_count.txt\", header = TRUE, sep = \"\\t\")\n", outdir);
   fprintf(gfp1,"r1_rs <- read.table(\"%s/Read1_rs_SN_MM_count.txt\", header = TRUE, sep = \"\\t\")\n", outdir);
   fprintf(gfp1,"r2_fs <- read.table(\"%s/Read2_fs_SN_MM_count.txt\", header = TRUE, sep = \"\\t\")\n", outdir);
   fprintf(gfp1,"r2_rs <- read.table(\"%s/Read2_rs_SN_MM_count.txt\", header = TRUE, sep = \"\\t\")\n", outdir);
   
   fprintf(gfp1,"pval <- paste(\"chi-square Read1_fs vs Read2_fs =\" ,chisq.test(cbind(r1_fs$no_of_reads, r2_fs$no_of_reads))$p.value)\n");
   fprintf(gfp1,"write(pval, file=\"%s/statistical_test.txt\", append=TRUE)\n", outdir);
   fprintf(gfp1,"pval <- paste(\"chi-square Read1_rs vs Read2_rs =\" ,chisq.test(cbind(r1_rs$no_of_reads, r2_rs$no_of_reads))$p.value)\n");
   fprintf(gfp1,"write(pval, file=\"%s/statistical_test.txt\", append=TRUE)\n", outdir);
   
   fprintf(gfp1,"pval <- paste(\"Kruskal-Wallis Read1_fs vs Read2_fs =\" ,kruskal.test(r1_fs$no_of_reads, r2_fs$no_of_reads)$p.value)\n");
   fprintf(gfp1,"write(pval, file=\"%s/statistical_test.txt\", append=TRUE)\n", outdir);
   fprintf(gfp1,"pval <- paste(\"Kruskal-Wallis Read1_rs vs Read2_rs =\" ,kruskal.test(r1_rs$no_of_reads, r2_rs$no_of_reads)$p.value)\n");
   fprintf(gfp1,"write(pval, file=\"%s/statistical_test.txt\", append=TRUE)\n", outdir);
   ////*updated on 19Apr2021 end*///
   
   fprintf(gfp1,"isizeplot <- read.table(\"%s/Overall_Insertsize_dist.txt\",header = TRUE, sep = \"\\t\")\n", outdir);
   fprintf(gfp1,"png(file=\"%s/Insertsize.png\")\n", outdir);
   fprintf(gfp1,"ggplot(data = isizeplot, aes(x = Isize, y = No_of_reads)) +  geom_bar(stat=\"identity\", fill=\"deeppink3\", width=0.3) + xlab(\"insert size (bp)\") + ylab(\"no of reads\")\n");
   fprintf(gfp1,"dev.off()\n");
   fprintf(gfp1,"gcplot <- read.table(\"%s/Overall_GC_dist.txt\",header = TRUE, sep = \"\\t\")\n", outdir);
   fprintf(gfp1,"png(file=\"%s/GC.png\")\n", outdir);
   fprintf(gfp1,"ggplot(data = gcplot, aes(x = GCpct, y = pct_of_reads)) +  geom_line(color=\"royalblue\") + expand_limits(x = 0, y = 0) + geom_point(color=\"gray30\",size =1) + xlab(\"GC content (%%) \") + ylab(\"%% of reads\")\n");
   fprintf(gfp1,"dev.off()\n");
   fprintf(gfp1,"mapqplot <- read.table(\"%s/Overall_Mappingquality_bin_dist.txt\",header = TRUE, sep = \"\t\")\n", outdir);
   fprintf(gfp1,"png(file=\"%s/Mapping_quality.png\")\n", outdir);
   fprintf(gfp1,"ggplot(data = mapqplot, aes(x = 2, y = Pct_reads, fill = MapQ)) +  geom_bar(stat=\"identity\", color=\"white\") + coord_polar(theta = \"y\") + theme_void() + xlim(c(0.5, 2.5)) + scale_fill_manual(values = c(\"red2\",\"skyblue\",\"skyblue2\",\"deepskyblue\",\"deepskyblue3\",\"dodgerblue3\",\"dodgerblue4\")) + scale_y_reverse()\n");
   fprintf(gfp1,"dev.off()\n");
   fprintf(gfp1,"baseqplot <- read.table(\"%s/Overall_Basequality_dist.txt\",header = TRUE, sep = \"\\t\")\n", outdir);
   fprintf(gfp1,"png(file=\"%s/Basequality.png\")\n", outdir);
   fprintf(gfp1,"ggplot(data = baseqplot, aes(x = BaseQ, y = Basecounts)) +  geom_bar(stat=\"identity\", fill=\"cyan3\", width = 0.75) + xlab(\"base quality\") + ylab(\"base counts\")\n");
   fprintf(gfp1,"dev.off()\n");
   fprintf(gfp1,"percybaseq <- read.table(\"%s/Overall_Percycle_avg_baseQ.txt\",header = TRUE, sep = \"\\t\")\n", outdir);
   fprintf(gfp1,"png(file=\"%s/Percycle_qv.png\")\n", outdir);
   fprintf(gfp1,"ggplot(data = percybaseq, aes(x = Readpos, y = Avg_baseQ)) +  geom_line(linetype=\"dotted\") + expand_limits(x = 0, y = 0) + geom_point(color=\"red\") + xlab(\"read position (bp)\") + ylab(\"quality (average)\")\n");
   fprintf(gfp1,"dev.off()\n");
   fprintf(gfp1,"percyneclu <- read.table(\"%s/Overall_Perbase_seq-content.txt\",header = TRUE, sep = \"\\t\")\n", outdir);
   fprintf(gfp1,"png(file=\"%s/Perbase-seq-content.png\")\n", outdir);
   fprintf(gfp1,"ggplot(data = percyneclu, aes(x = Readpos, y = Pct, group = Base)) +  geom_line(aes(color=Base)) + scale_color_manual(values=c(\"gray30\",\"darkorange3\",\"springgreen4\",\"red1\",\"deepskyblue3\")) + geom_point(aes(color=Base),size=0.3) + xlab(\"read position (bp)\") + ylab(\"nucleotide content (%%) \") + expand_limits(y=c(0,100))\n");
   fprintf(gfp1,"dev.off()\n");
   fclose(gfp1);
   *outfn=0;
   sprintf(outfn, "R CMD BATCH %s/plotsmake.r", outdir);
   system(outfn);
   //printf("Command :: %s\n", outfn);
   //printf("Plots Created\n");
   *outfn=0;
   sprintf(outfn, "rm -f %s/plotsmake.r", outdir);
   system(outfn);
   system("rm -f plotsmake.*");
   free(outfn);
 }

 void makeplots_v1(char *outdir, char *prefix)
 {
   FILE *gfp1;
   char *outfn;
   outfn = calloc(strlen(outdir) + 500, 1);
   //printf("%s	%s\n",outdir, prefix);
   sprintf(outfn, "%s/makeplots.r", outdir);
   gfp1=fopen(outfn,"w");
   *outfn=0;
   sprintf(outfn, "%s/%s", outdir, prefix);
   //printf("Output file with path = %s\n", outfn);
   fprintf(gfp1,"library (\"ggplot2\")\n");

   fprintf(gfp1,"library (\"gridExtra\")\n");

   fprintf(gfp1,"bplot <- read.table(\"%s_ReadposTotal_SN_MM.txt\",header = TRUE, sep = \"\t\")\n", outfn);
   fprintf(gfp1,"qvplot <- read.table(\"%s_Readpos_MMQV.txt\", header = TRUE, sep = \"\t\")\n", outfn);
   fprintf(gfp1,"hmap <- read.table(\"%s_ReadposBaseChange_SN_MM.txt\",header = TRUE, sep = \"\t\")\n", outfn);
   fprintf(gfp1,"hmap_order <- c(\"A>C\",\"T>G\",\"A>G\",\"T>C\",\"A>T\",\"T>A\",\"C>A\",\"G>T\",\"C>G\",\"G>C\",\"C>T\",\"G>A\")\n");

   fprintf(gfp1,"g1 <- ggplot(data = bplot, aes(x = Pos, y = pct_of_reads)) +  geom_bar(stat=\"identity\", colour = \"gray95\", fill=\"deepskyblue3\") + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin=unit(c(0.1,0.1,-0.15,0.1), \"cm\")) + ylab(\"reads having\n mismatches (%%)\") + labs(title=\"%s\") + theme(plot.title = element_text(hjust = 1))\n",prefix);

   fprintf(gfp1,"g2 <- ggplot(data = qvplot, aes(x = Pos, y = Avg_qv)) +  geom_line(linetype=\"dotted\") + geom_point(color=\"violetred1\") + ylab(\"missmatch base\nquality (average)\") + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.margin=unit(c(0.05,0.18,0.1,0.2), \"cm\"))\n");

   fprintf(gfp1,"g3 <- ggplot(data = hmap, aes(x = Pos, y = Change)) +  geom_tile(aes(fill = Value), colour = \"white\",size =0.3) + scale_fill_gradient(\"Change (%%)\",low = \"white\",high = \"red2\") + theme(legend.position = \"bottom\", legend.direction = \"horizontal\",plot.margin=unit(c(-0.15,0.1,0.1,0.1), \"cm\")) + xlab(\"Read position (bp)\") + ylab(\"change\") + scale_y_discrete(limits = hmap_order)\n");

   fprintf(gfp1,"gg1 <- ggplot_gtable(ggplot_build(g1))\n");
   fprintf(gfp1,"gg2 <- ggplot_gtable(ggplot_build(g2))\n");
   fprintf(gfp1,"gg3 <- ggplot_gtable(ggplot_build(g3))\n");

   fprintf(gfp1,"maxWidth = grid::unit.pmax(gg1$widths[2:5], gg2$widths[2:5], gg3$widths[2:5])\n");
   fprintf(gfp1,"gg1$widths[2:5] <- as.list(maxWidth)\n");
   fprintf(gfp1,"gg2$widths[2:5] <- as.list(maxWidth)\n");
   fprintf(gfp1,"gg3$widths[2:5] <- as.list(maxWidth)\n");
   fprintf(gfp1,"png(file=\"%s_Com_plots_U.png\",width=10, height=7, units=\"in\", res=100)\n", outfn);
   fprintf(gfp1,"grid.arrange(gg1, gg2, gg3)\n");
   fprintf(gfp1,"dev.off()\n");

   fprintf(gfp1,"mmbplot <- read.table(\"%s_SN_MM.txt\",header = TRUE, sep = \"\t\")\n", outfn);
   fprintf(gfp1,"mmodr <- c(\"1\",\"2\",\"3\",\"4\",\">=5\")\n");
   fprintf(gfp1,"png(file=\"%s_Mismatchplot_U.png\",res=100)\n", outfn);
   fprintf(gfp1,"ggplot(data = mmbplot, aes(x = Missmatch, y = pct_of_reads)) +  geom_bar(stat=\"identity\", fill=c(\"#FF9999\",\"#FF6666\",\"#FF3333\",\"#FF0000\",\"#CC0033\")) + ylab(\"%% of reads\") + xlab(\"mismatch count\") + scale_x_discrete(limits = mmodr) + labs(title=\"%s\") + theme(plot.title = element_text(hjust = 1))\n",prefix);
   fprintf(gfp1,"dev.off()\n");

   fprintf(gfp1,"bcng <- read.table(\"%s_BaseChange_SN_MM.txt\",header = TRUE, sep = \"\t\")\n", outfn);
   fprintf(gfp1,"bcng_order <- c(\"A>C\",\"A>G\",\"A>T\",\"C>A\",\"C>G\",\"C>T\",\"G>A\",\"G>C\",\"G>T\",\"T>A\",\"T>C\",\"T>G\")\n");

   fprintf(gfp1,"bcngQV <- read.table(\"%s_BaseChangeQV_SN_MM.txt\",header = TRUE, sep = \"\t\")\n", outfn);
   fprintf(gfp1,"bcngQV_order <- c(\"A>C\",\"A>G\",\"A>T\",\"C>A\",\"C>G\",\"C>T\",\"G>A\",\"G>C\",\"G>T\",\"T>A\",\"T>C\",\"T>G\")\n");
  
   fprintf(gfp1,"mmflnkshmap <- read.table(\"%s_BaseChange_SN_MM_flanks.txt\",header = TRUE, sep = \"\t\")\n", outfn);
   
   fprintf(gfp1,"g1 <- ggplot(data = bcng, aes(x = Change, y = pct_of_bases)) +  geom_bar(stat=\"identity\", fill=\"darkorange\") + ylab(\"change (%%)\") + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin=unit(c(0.1,2.9,-0.15,0.1), \"cm\")) + scale_x_discrete(limits = bcng_order) + labs(title=\"%s\") + theme(plot.title = element_text(hjust = 1))\n",prefix);

   fprintf(gfp1,"g2 <- ggplot(data = bcngQV, aes(x = Change, y = QV,fill=factor(Range, levels=c(\">30\",\"21-30\",\"11-20\",\"<=10\")))) +  geom_bar(stat=\"identity\") + scale_fill_manual(\"Base quality\",values = c(\"green4\",\"green3\",\"seagreen2\",\"red1\")) + ylab(\"basequality bins (%%)\") + labs(\"QV\") + scale_x_discrete(limits = bcngQV_order) + theme(plot.margin=unit(c(-0.15,0.01,-0.15,0.1), \"cm\")) + theme(axis.title.x=element_blank(),axis.text.x=element_blank())\n");

   fprintf(gfp1,"g3 <- ggplot(data = mmflnkshmap, aes(y = Flanks, x = Change)) +  geom_tile(aes(fill = Pct), colour = \"white\",size =0.3) + scale_fill_gradient(\"flanks  (%%)\",low = \"gray95\",high = \"deeppink3\") + xlab(\"change\") + ylab(\"flanking bases (5'_3')\") + scale_x_discrete(limits = bcng_order) + theme(plot.margin=unit(c(-0.05,0.28,0.01,0.12), \"cm\"))\n");
   
   fprintf(gfp1,"gg1 <- ggplot_gtable(ggplot_build(g1))\n");
   fprintf(gfp1,"gg2 <- ggplot_gtable(ggplot_build(g2))\n");
   fprintf(gfp1,"gg3 <- ggplot_gtable(ggplot_build(g3))\n");
   fprintf(gfp1,"maxWidth = grid::unit.pmax(gg1$widths[2:5], gg2$widths[2:5], gg3$widths[2:5])\n");
   fprintf(gfp1,"gg1$widths[2:5] <- as.list(maxWidth)\n");
   fprintf(gfp1,"gg2$widths[2:5] <- as.list(maxWidth)\n");
   fprintf(gfp1,"gg3$widths[2:5] <- as.list(maxWidth)\n");
   fprintf(gfp1,"png(file=\"%s_Basechange_and_Quality_U.png\")\n", outfn);
   fprintf(gfp1,"grid.arrange(gg1, gg2, gg3)\n");
   fprintf(gfp1,"dev.off()\n");

   fprintf(gfp1,"clip <- read.table(\"%s_CLIP_MM.txt\",header = TRUE, sep = \"\t\")\n", outfn);
   fprintf(gfp1,"clip_order <- c(\"1-5\",\"6-10\",\"11-15\",\"16-20\",\">20\")\n");

   fprintf(gfp1,"png(file=\"%s_Clipped_plot_U.png\",res=100)\n", outfn);
   fprintf(gfp1,"ggplot(data = clip , aes(x = Bases, y = pct_of_reads, fill=factor(End,levels=c(\"5prime\",\"3prime\")))) +  geom_bar(stat=\"identity\",position = \"dodge\") + scale_x_discrete(limits = clip_order) + scale_fill_manual(values = c(\"gray30\",\"brown2\")) + theme(legend.title = element_blank()) + xlab(\"Clipped bases\") + ylab(\"%% of reads\") + labs(title=\"%s\") + theme(plot.title = element_text(hjust = 1))\n",prefix);
   fprintf(gfp1,"dev.off()\n");

      
   fprintf(gfp1,"indel <- read.table(\"%s_INDEL_MM.txt\",header = TRUE, sep = \"\t\")\n", outfn);
   fprintf(gfp1,"indel_order <- c(\">=-5\",\"-4\",\"-3\",\"-2\",\"-1\",\"0\",\"+1\",\"+2\",\"+3\",\"+4\",\">=+5\")\n");
   //fprintf(gfp1,"flnkshmap_indel <- read.table(\"%s_INDEL_MM_flanks.txt\",header = TRUE, sep = \"\t\")\n", outfn);
   
   
   fprintf(gfp1,"png(file=\"%s_Indel_plot_U.png\")\n", outfn);
   fprintf(gfp1,"ggplot(data = indel, aes(x = factor(Len), y = pct_of_reads, fill=factor(Type, levels=c(\"C\",\"G\",\"T\",\"A\",\"Multi-nucleotide\")))) +  geom_bar(stat=\"identity\") + scale_fill_manual(values = c(\"lightsalmon4\",\"lightseagreen\",\"brown2\",\"gray30\",\"darkorange\")) + theme(legend.title = element_blank()) + scale_x_discrete(limits = indel_order) + ylab(\"%% of reads\") + labs(title=\"%s\") + xlab(\"Deletion                                                    Insertion\") + theme(plot.title = element_text(hjust = 1))\n",prefix);

   //fprintf(gfp1,"ggplot(data = indel, aes(x = factor(Len), y = pct_of_reads, fill=factor(Type, levels=c(\"C\",\"G\",\"T\",\"A\",\"Multi-nucleotide\")))) +  geom_bar(stat=\"identity\") + xlab(\"Deletion                                                    Insertion\") + scale_fill_manual(values = c(\"lightsalmon4\",\"lightseagreen\",\"brown2\",\"gray30\",\"darkorange\")) + theme(legend.position = \"bottom\", legend.direction = \"horizontal\") + theme(legend.title = element_blank()) + scale_x_discrete(limits = indel_order) + ylab(\"%% of reads\")\n");
   
   //fprintf(gfp1,"g2 <- ggplot(data = flnkshmap_indel, aes(y = Flanks, x = Length)) +  geom_tile(aes(fill = Pct), colour = \"white\",size =0.3) + scale_fill_gradient(\"flanks  (%%)\",low = \"gray95\",high = \"deeppink3\") + xlab(\"Deletion                                                    Insertion\") + ylab(\"flanking bases (5'_3')\") + scale_x_discrete(limits = indel_order) + theme(plot.margin=unit(c(-0.40,1.30,0.2,0.15), \"cm\"))\n");
   
   //fprintf(gfp1,"gg1 <- ggplot_gtable(ggplot_build(g1))\n");
   //fprintf(gfp1,"gg2 <- ggplot_gtable(ggplot_build(g2))\n");
   //fprintf(gfp1,"maxWidth = grid::unit.pmax(gg1$widths[2:5], gg2$widths[2:5])\n");
   //fprintf(gfp1,"gg1$widths[2:5] <- as.list(maxWidth)\n");
   //fprintf(gfp1,"gg2$widths[2:5] <- as.list(maxWidth)\n");
   //fprintf(gfp1,"png(file=\"%s_Indel_plot_U.png\")\n", outfn);
   //fprintf(gfp1,"grid.arrange(gg1, gg2);\n");
   fprintf(gfp1,"dev.off()\n");
   fclose(gfp1);
   *outfn=0;
   sprintf(outfn, "R CMD BATCH %s/makeplots.r", outdir);
   system(outfn);
   *outfn=0;
   sprintf(outfn, "rm -f %s/makeplots.r", outdir);
   system(outfn);
   system("rm -f makeplots.*");
   system("rm -f Rplots.pdf");
   
   free(outfn);
 }

void crthtml(paramtr *prmtr, char *tme, log_v *OA, log_v *R1, log_v *R2)
{

  FILE *fp;
  int i;
  float misize, mbasq, mmapq, mrdln;
  char ch, *outfn=0;
  //printf("I am in Html\n");
  outfn = calloc(strlen(prmtr->outdir) + 500, 1);
  sprintf(outfn, "%s/Overall_mapping_summary.log", prmtr->outdir);
  //printf("Log file = %s\n",outfn);
  fp=fopen(outfn,"r");
  for(i=1;i<=61;i++) { while((ch=fgetc(fp))!='\n'); }
  while((ch=fgetc(fp))!='=');
  ch=fgetc(fp);
  fscanf(fp,"%f",&misize);
  ch=fgetc(fp);
  
  while((ch=fgetc(fp))!='=');
  ch=fgetc(fp);
  fscanf(fp,"%f",&mbasq);
  ch=fgetc(fp);
  
  while((ch=fgetc(fp))!='=');
  ch=fgetc(fp);
  fscanf(fp,"%f",&mmapq);
  ch=fgetc(fp);
  
  while((ch=fgetc(fp))!='=');
  ch=fgetc(fp);
  fscanf(fp,"%f",&mrdln);
  ch=fgetc(fp);
  
  fclose(fp);
  *outfn = 0;
  sprintf(outfn, "%s/Bamqc.html", prmtr->outdir);
  fp=fopen(outfn,"w");
	fprintf(fp,"<!DOCTYPE html>\n");
	fprintf(fp,"<html lang=\"en\">\n");
	fprintf(fp,"<head>\n");
	fprintf(fp,"<meta charset=\"utf-8\">\n");
	fprintf(fp,"<title>Report(bamqc)</title>\n");
	fprintf(fp,"<style>\n");
	fprintf(fp,"    body{        \n");
	fprintf(fp,"        padding-top: 70px;\n");
	fprintf(fp,"        padding-bottom: 40px;\n");
	fprintf(fp,"        padding-left: 20em;\n");
	fprintf(fp,"        \n");
	fprintf(fp,"    }\n");
	fprintf(fp,"    .container{\n");
	fprintf(fp,"        width: 90%%;\n");
	fprintf(fp,"        margin: 0 auto; /* Center the DIV horizontally */\n");
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
//	fprintf(fp,"        <div class=\"container\">Local QC</div>\n");
        fprintf(fp,"<div class=\"headertitle\">\n");
        fprintf(fp,"<span style=\"color:#FF5733; font-size:125%%; float: left; padding-left: 20px;\"><b>Map</b></span>\n");
        fprintf(fp,"<span style=\"color:#16A085; font-size:125%%; float: left;\"><b>Insights</b></span>\n");
        fprintf(fp,"<mid style=\"color:white; font-size:100%%; float: right; padding-right: 20px;\">Mapinsights bamqc report</mid>\n");
        fprintf(fp,"</div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<ul class=\"navbar\">\n");
	fprintf(fp,"<h2 style=\"color:#797D7F;\">CONTENTS</h2>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Parameters & inputs\" style=\"color:#1C2833;\" >Commands and parameters</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Summary statistics\" style=\"color:#1C2833;\" >Summary statistics</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Base quality distribution\" style=\"color:#1C2833;\" >Base quality distribution</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Per base average quality\" style=\"color:#1C2833;\" >Per base average quality</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Per base sequence content\" style=\"color:#1C2833;\" >Per base sequence content</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#GC content distribution\" style=\"color:#1C2833;\" >GC content distribution</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Mapping quality distribution\" style=\"color:#1C2833;\" >Mapping quality profile</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Insert size histogram\" style=\"color:#1C2833;\" >Insert size histogram</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Mismatch counts\" style=\"color:#1C2833;\" >Mismatch counts</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Basechange & quality\" style=\"color:#1C2833;\" >Basechange & quality</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Per cycle mismatchs & quality status\" style=\"color:#1C2833;\" >Per cycle mismatchs & quality status</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Deletion and Insertion status\" style=\"color:#1C2833;\" >Deletion and Insertion status</a></li>\n");
	fprintf(fp,"<li class=\"toctree-l1\" onmouseover=\"style.fontWeight = 'bold'\" onmouseout=\"style.fontWeight = 'normal'\"><a class=\"reference internal\" href=\"#Read clipping status\" style=\"color:#1C2833;\" >Read clipping status</a></li>\n");
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
	fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Command line</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>mapinsights %s</td>\n", prmtr->comnd);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Parameters</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Alignment file</td>\n");
	fprintf(fp,"<td class=column2>%s</td>\n", prmtr->bam_file);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Reference file</td>\n");
	fprintf(fp,"<td class=column2>%s</td>\n", prmtr->ref_file);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Bed file</td>\n");
	if(prmtr->bed_file) {fprintf(fp,"<td class=column2>%s</td>\n", prmtr->bed_file);}
	else {fprintf(fp,"<td class=column2>Not specified</td>\n");}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	/*fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Threads</td>\n");
	fprintf(fp,"<td class=column2>%d</td>\n", prmtr->thrd);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");*/
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Other information</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"\n");
        fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>No of contigs in alignment file</td>\n");
	fprintf(fp,"<td class=column2>%d</td>\n", prmtr->no_of_contigs);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
        fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>No of @RG tag</td>\n");
        if(prmtr->no_of_rgs) {fprintf(fp,"<td class=column2>%d</td>\n", prmtr->no_of_rgs);}
        else {fprintf(fp,"<td class=column2>Not specified</td>\n");}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Sample name [@RG SM]</td>\n");
        if(prmtr->sample_name){fprintf(fp,"<td class=column2>%s</td>\n", prmtr->sample_name);}
        else {fprintf(fp,"<td class=column2>Not specified</td>\n");}
        fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
        
        fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Library [@RG LB]</td>\n");
        if(prmtr->library) {fprintf(fp,"<td class=column2>%s</td>\n", prmtr->library);}
        else {fprintf(fp,"<td class=column2>Not specified</td>\n");}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Identifier [@RG ID]</td>\n");
        if(prmtr->no_of_rgs) 
        {
           fprintf(fp,"<td class=column2>");
           for(i=0;i<prmtr->no_of_rgs;i++) {fprintf(fp,"%s ", prmtr->sId[i]);}
           fprintf(fp,"</td>\n");
        }
        else {fprintf(fp,"<td class=column2>Not specified</td>\n");}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	
	if(prmtr->xrg_file && prmtr->xrg_cnt > 0)
	{
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Exclude identifier list</td>\n");
        if(prmtr->no_of_rgs) 
        {
           fprintf(fp,"<td class=column2>");
           for(i=0;i<prmtr->xrg_cnt;i++) {fprintf(fp,"%s ", prmtr->x_rg[i]);}
           fprintf(fp,"</td>\n");
        }
        else {fprintf(fp,"<td class=column2>Not specified</td>\n");}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	}
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Analysis date</td>\n");
	fprintf(fp,"<td class=column2>%s</td>\n", tme);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Homopolymer size</td>\n");
	fprintf(fp,"<td class=column2>5</td>\n");
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Reference size (bp)</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", prmtr->ref_size);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Targeted region size (bp)</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", prmtr->trgtsize);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,"<h2 id=\"Summary statistics\" style=\"background-color:#E5E7E9;\"><u>Summary statistics</u></h2>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Global counts</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>No of reads in input file</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->glog.all_rd_cnt);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Unmapped reads</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.unmap, ((float)OA->glog.unmap/(float)OA->glog.all_rd_cnt)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	if(prmtr->trgtsize > 0){
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>No of reads (out of targets)</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", (OA->glog.all_rd_cnt - (OA->glog.tot_rd_cnt + OA->glog.unmap)));
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	}
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=mapsummary>\n");
	if(prmtr->trgtsize==0) {fprintf(fp,"<h3>Mapping logs (Globals)</h3>\n");}
	else {fprintf(fp,"<h3>Mapping logs (Inside of targeted regions)</h3>\n");}
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1></td>\n");
	fprintf(fp,"<td class=column2> <b>Overall</b></td>\n");
	fprintf(fp,"<td class=column3> <b>Read1</b></td>\n");
	fprintf(fp,"<td class=column4> <b>Read2</b></td>\n");
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mapped-reads</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->glog.tot_rd_cnt); 
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.tot_rd_cnt, ((float)R1->glog.tot_rd_cnt/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.tot_rd_cnt, ((float)R2->glog.tot_rd_cnt/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mapped-read1</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.r1_cnt, ((float)OA->glog.r1_cnt/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.r1_cnt, ((float)R1->glog.r1_cnt/(float)OA->glog.r1_cnt)*100);
	fprintf(fp,"<td class=column4>-</td>\n");
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mapped-read2</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.r2_cnt, ((float)OA->glog.r2_cnt/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"<td class=column3>-</td>\n");
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.r2_cnt, ((float)R2->glog.r2_cnt/(float)OA->glog.r2_cnt)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mapped-forward</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", (OA->glog.tot_rd_cnt - OA->glog.rev_cnt), ((float)(OA->glog.tot_rd_cnt - OA->glog.rev_cnt)/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", (R1->glog.tot_rd_cnt - R1->glog.rev_cnt), ((float)(R1->glog.tot_rd_cnt - R1->glog.rev_cnt)/(float)(OA->glog.tot_rd_cnt - OA->glog.rev_cnt))*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", (R2->glog.tot_rd_cnt - R2->glog.rev_cnt), ((float)(R2->glog.tot_rd_cnt - R2->glog.rev_cnt)/(float)(OA->glog.tot_rd_cnt - OA->glog.rev_cnt))*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mapped-reverse</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.rev_cnt, ((float)OA->glog.rev_cnt/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.rev_cnt, ((float)R1->glog.rev_cnt/(float)OA->glog.rev_cnt)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.rev_cnt, ((float)R2->glog.rev_cnt/(float)OA->glog.rev_cnt)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mapped-pair</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.map_pr, ((float)OA->glog.map_pr/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.map_pr, ((float)R1->glog.map_pr/(float)OA->glog.map_pr)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.map_pr, ((float)R2->glog.map_pr/(float)OA->glog.map_pr)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mapped-properpair</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.propr_pr, ((float)OA->glog.propr_pr/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.propr_pr, ((float)R1->glog.propr_pr/(float)OA->glog.propr_pr)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.propr_pr, ((float)R2->glog.propr_pr/(float)OA->glog.propr_pr)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Secondary-alignments</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.scndry_cnt, ((float)OA->glog.scndry_cnt/(float)OA->glog.tot_rd_cnt)*100);
	if(OA->glog.scndry_cnt == 0) {
	fprintf(fp,"<td class=column3>0 (0.00%%)</td>\n");
	fprintf(fp,"<td class=column4>0 (0.00%%)</td>\n");
	} else {
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.scndry_cnt, ((float)R1->glog.scndry_cnt/(float)OA->glog.scndry_cnt)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.scndry_cnt, ((float)R2->glog.scndry_cnt/(float)OA->glog.scndry_cnt)*100);}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Supplementary-alignments</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.supply_cnt, ((float)OA->glog.supply_cnt/(float)OA->glog.tot_rd_cnt)*100);
	if(OA->glog.supply_cnt == 0) {
	fprintf(fp,"<td class=column3>0 (0.00%%)</td>\n");
	fprintf(fp,"<td class=column4>0 (0.00%%)</td>\n");
	} else {
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.supply_cnt, ((float)R1->glog.supply_cnt/(float)OA->glog.supply_cnt)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.supply_cnt, ((float)R2->glog.supply_cnt/(float)OA->glog.supply_cnt)*100);}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>QC-failed</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.qc_fail, ((float)OA->glog.qc_fail/(float)OA->glog.tot_rd_cnt)*100);
	if(OA->glog.qc_fail == 0) {
	fprintf(fp,"<td class=column3>0 (0.00%%)</td>\n");
	fprintf(fp,"<td class=column4>0 (0.00%%)</td>\n");
	} else {
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.qc_fail, ((float)R1->glog.qc_fail/(float)OA->glog.qc_fail)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.qc_fail, ((float)R2->glog.qc_fail/(float)OA->glog.qc_fail)*100);}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Strand-ratio (F:R)</td>\n");
	fprintf(fp,"<td class=column2>%.2f:%.2f</td>\n", (float)(OA->glog.tot_rd_cnt - OA->glog.rev_cnt)/(float)OA->glog.tot_rd_cnt,(float)OA->glog.rev_cnt/(float)OA->glog.tot_rd_cnt);
	fprintf(fp,"<td class=column3>%.2f:%.2f</td>\n", (float)(R1->glog.tot_rd_cnt - R1->glog.rev_cnt)/(float)R1->glog.tot_rd_cnt,(float)R1->glog.rev_cnt/(float)R1->glog.tot_rd_cnt);
	fprintf(fp,"<td class=column4>%.2f:%.2f</td>\n", (float)(R2->glog.tot_rd_cnt - R2->glog.rev_cnt)/(float)R2->glog.tot_rd_cnt,(float)R2->glog.rev_cnt/(float)R2->glog.tot_rd_cnt);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Clipped summary</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Softclip-events</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->midc.sftclp);
	if(OA->midc.sftclp == 0) {
	fprintf(fp,"<td class=column3>0 (0.00%%)</td>\n");
	fprintf(fp,"<td class=column4>0 (0.00%%)</td>\n");
	} else {
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->midc.sftclp, ((float)R1->midc.sftclp/(float)OA->midc.sftclp)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->midc.sftclp, ((float)R2->midc.sftclp/(float)OA->midc.sftclp)*100);}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Softclipped basecounts</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->midc.sftclp_bs);
	if(OA->midc.sftclp_bs == 0) {
	fprintf(fp,"<td class=column3>0 (0.00%%)</td>\n");
	fprintf(fp,"<td class=column4>0 (0.00%%)</td>\n");
	} else {
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->midc.sftclp_bs, ((float)R1->midc.sftclp_bs/(float)OA->midc.sftclp_bs)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->midc.sftclp_bs, ((float)R2->midc.sftclp_bs/(float)OA->midc.sftclp_bs)*100);}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Hardclip-events</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->midc.hrdclp);
	if(OA->midc.hrdclp == 0) {
	fprintf(fp,"<td class=column3>0 (0.00%%)</td>\n");
	fprintf(fp,"<td class=column4>0 (0.00%%)</td>\n");
	} else {
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->midc.hrdclp, ((float)R1->midc.hrdclp/(float)OA->midc.hrdclp)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->midc.hrdclp, ((float)R2->midc.hrdclp/(float)OA->midc.hrdclp)*100);}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Hardclipped basecounts</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->midc.hrdclp_bs);
	if(OA->midc.hrdclp_bs == 0) {
	fprintf(fp,"<td class=column3>0 (0.00%%)</td>\n");
	fprintf(fp,"<td class=column4>0 (0.00%%)</td>\n");
	} else {
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->midc.hrdclp_bs, ((float)R1->midc.hrdclp_bs/(float)OA->midc.hrdclp_bs)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->midc.hrdclp_bs, ((float)R2->midc.hrdclp_bs/(float)OA->midc.hrdclp_bs)*100);}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Adapter summary</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Illumina-Adapter</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->adptr.IADptr);
	if(OA->adptr.IADptr == 0) {
	fprintf(fp,"<td class=column3>0 (0.00%%)</td>\n");
	fprintf(fp,"<td class=column4>0 (0.00%%)</td>\n");
	} else {
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->adptr.IADptr, ((float)R1->adptr.IADptr/(float)OA->adptr.IADptr)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->adptr.IADptr, ((float)R2->adptr.IADptr/(float)OA->adptr.IADptr)*100);}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Illumina-PCR-primer</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->adptr.IPCRprm);
	if(OA->adptr.IPCRprm == 0) {
	fprintf(fp,"<td class=column3>0 (0.00%%)</td>\n");
	fprintf(fp,"<td class=column4>0 (0.00%%)</td>\n");
	} else {
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->adptr.IPCRprm, ((float)R1->adptr.IPCRprm/(float)OA->adptr.IPCRprm)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->adptr.IPCRprm, ((float)R2->adptr.IPCRprm/(float)OA->adptr.IPCRprm)*100);}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Nextera-Transposase-Sequence</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->adptr.NTS);
	if(OA->adptr.NTS == 0) {
	fprintf(fp,"<td class=column3>0 (0.00%%)</td>\n");
	fprintf(fp,"<td class=column4>0 (0.00%%)</td>\n");
	} else {
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->adptr.NTS, ((float)R1->adptr.NTS/(float)OA->adptr.NTS)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->adptr.NTS, ((float)R2->adptr.NTS/(float)OA->adptr.NTS)*100);}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Base content</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>A</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.AAcnt, ((double)(OA->glog.AAcnt)/(double)(OA->glog.AAcnt+OA->glog.TTcnt+OA->glog.GGcnt+OA->glog.CCcnt+OA->glog.NNcnt)*100));
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.AAcnt, ((float)R1->glog.AAcnt/(float)OA->glog.AAcnt)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.AAcnt, ((float)R2->glog.AAcnt/(float)OA->glog.AAcnt)*100);
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>T</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.TTcnt, ((double)(OA->glog.TTcnt)/(double)(OA->glog.AAcnt+OA->glog.TTcnt+OA->glog.GGcnt+OA->glog.CCcnt+OA->glog.NNcnt)*100));
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.TTcnt, ((float)R1->glog.TTcnt/(float)OA->glog.TTcnt)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.TTcnt, ((float)R2->glog.TTcnt/(float)OA->glog.TTcnt)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>G</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.GGcnt, ((double)(OA->glog.GGcnt)/(double)(OA->glog.AAcnt+OA->glog.TTcnt+OA->glog.GGcnt+OA->glog.CCcnt+OA->glog.NNcnt)*100));
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.GGcnt, ((float)R1->glog.GGcnt/(float)OA->glog.GGcnt)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.GGcnt, ((float)R2->glog.GGcnt/(float)OA->glog.GGcnt)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>C</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.CCcnt, ((double)(OA->glog.CCcnt)/(double)(OA->glog.AAcnt+OA->glog.TTcnt+OA->glog.GGcnt+OA->glog.CCcnt+OA->glog.NNcnt)*100));
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.CCcnt, ((float)R1->glog.CCcnt/(float)OA->glog.CCcnt)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.CCcnt, ((float)R2->glog.CCcnt/(float)OA->glog.CCcnt)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>N</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.NNcnt, ((double)(OA->glog.NNcnt)/(double)(OA->glog.AAcnt+OA->glog.TTcnt+OA->glog.GGcnt+OA->glog.CCcnt+OA->glog.NNcnt)*100));
	if(OA->glog.NNcnt == 0) {
	fprintf(fp,"<td class=column3>0 (0.00%%)</td>\n");
	fprintf(fp,"<td class=column4>0 (0.00%%)</td>\n");
	} else {
	fprintf(fp,"<td class=column3>%ld (%.2f%%)</td>\n", R1->glog.NNcnt, ((float)R1->glog.NNcnt/(float)OA->glog.NNcnt)*100);
	fprintf(fp,"<td class=column4>%ld (%.2f%%)</td>\n", R2->glog.NNcnt, ((float)R2->glog.NNcnt/(float)OA->glog.NNcnt)*100);}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>GC %%</td>\n");
	fprintf(fp,"<td class=column2>%.2f</td>\n", ((double)(OA->glog.GGcnt+OA->glog.CCcnt)/(double)(OA->glog.AAcnt+OA->glog.TTcnt+OA->glog.GGcnt+OA->glog.CCcnt)*100));
	fprintf(fp,"<td class=column3>%.2f</td>\n", ((double)(R1->glog.GGcnt+R1->glog.CCcnt)/(double)(R1->glog.AAcnt+R1->glog.TTcnt+R1->glog.GGcnt+R1->glog.CCcnt)*100));
	fprintf(fp,"<td class=column4>%.2f</td>\n", ((double)(R2->glog.GGcnt+R2->glog.CCcnt)/(double)(R2->glog.AAcnt+R2->glog.TTcnt+R2->glog.GGcnt+R2->glog.CCcnt)*100));
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Mismatches and indels</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Single nucleotide match</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", (long int)(OA->glog.AAcnt+OA->glog.TTcnt+OA->glog.GGcnt+OA->glog.CCcnt)-(OA->midc.basemmc));
	fprintf(fp,"<td class=column3>%ld</td>\n", (long int)(R1->glog.AAcnt+R1->glog.TTcnt+R1->glog.GGcnt+R1->glog.CCcnt)-(R1->midc.basemmc));
	fprintf(fp,"<td class=column4>%ld</td>\n", (long int)(R2->glog.AAcnt+R2->glog.TTcnt+R2->glog.GGcnt+R2->glog.CCcnt)-(R2->midc.basemmc));
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Single nucleotide mismatch</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->midc.basemmc);
	fprintf(fp,"<td class=column3>%ld</td>\n", R1->midc.basemmc);
	fprintf(fp,"<td class=column4>%ld</td>\n", R2->midc.basemmc);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Insertion-events</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->midc.ins_cnt);
	fprintf(fp,"<td class=column3>%ld</td>\n", R1->midc.ins_cnt);
	fprintf(fp,"<td class=column4>%ld</td>\n", R2->midc.ins_cnt);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Inserted basecount</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->midc.ins_cnt_bs);
	fprintf(fp,"<td class=column3>%ld</td>\n", R1->midc.ins_cnt_bs);
	fprintf(fp,"<td class=column4>%ld</td>\n", R2->midc.ins_cnt_bs);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Homopolymer Insertation</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->midc.polyIns);
	fprintf(fp,"<td class=column3>%ld</td>\n", R1->midc.polyIns);
	fprintf(fp,"<td class=column4>%ld</td>\n", R2->midc.polyIns);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Deletion-events</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->midc.del_cnt);
	fprintf(fp,"<td class=column3>%ld</td>\n", R1->midc.del_cnt);
	fprintf(fp,"<td class=column4>%ld</td>\n", R2->midc.del_cnt);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Deleted basecount</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->midc.del_cnt_bs);
	fprintf(fp,"<td class=column3>%ld</td>\n", R1->midc.del_cnt_bs);
	fprintf(fp,"<td class=column4>%ld</td>\n", R2->midc.del_cnt_bs);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Homopolymer Deletion</td>\n");
	fprintf(fp,"<td class=column2>%ld</td>\n", OA->midc.polyDel);
	fprintf(fp,"<td class=column3>%ld</td>\n", R1->midc.polyDel);
	fprintf(fp,"<td class=column4>%ld</td>\n", R2->midc.polyDel);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mismatch rate</td>\n");
	fprintf(fp,"<td class=column2>%.4f</td>\n", (float)(OA->midc.basemmc)/(float)(OA->glog.AAcnt+OA->glog.TTcnt+OA->glog.GGcnt+OA->glog.CCcnt));
	fprintf(fp,"<td class=column3>%.4f</td>\n", (float)(R1->midc.basemmc)/(float)(R1->glog.AAcnt+R1->glog.TTcnt+R1->glog.GGcnt+R1->glog.CCcnt));
	fprintf(fp,"<td class=column4>%.4f</td>\n", (float)(R2->midc.basemmc)/(float)(R2->glog.AAcnt+R2->glog.TTcnt+R2->glog.GGcnt+R2->glog.CCcnt));
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Insertion rate</td>\n");
	fprintf(fp,"<td class=column2>%.4f</td>\n", (float)(OA->midc.ins_cnt)/(float)(OA->glog.AAcnt+OA->glog.TTcnt+OA->glog.GGcnt+OA->glog.CCcnt));
	fprintf(fp,"<td class=column3>%.4f</td>\n", (float)(R1->midc.ins_cnt)/(float)(R1->glog.AAcnt+R1->glog.TTcnt+R1->glog.GGcnt+R1->glog.CCcnt));
	fprintf(fp,"<td class=column4>%.4f</td>\n", (float)(R2->midc.ins_cnt)/(float)(R2->glog.AAcnt+R2->glog.TTcnt+R2->glog.GGcnt+R2->glog.CCcnt));
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Deletion rate</td>\n");
	fprintf(fp,"<td class=column2>%.4f</td>\n", (float)(OA->midc.del_cnt)/(float)(OA->glog.AAcnt+OA->glog.TTcnt+OA->glog.GGcnt+OA->glog.CCcnt));
	fprintf(fp,"<td class=column3>%.4f</td>\n", (float)(R1->midc.del_cnt)/(float)(R1->glog.AAcnt+R1->glog.TTcnt+R1->glog.GGcnt+R1->glog.CCcnt));
	fprintf(fp,"<td class=column4>%.4f</td>\n", (float)(R2->midc.del_cnt)/(float)(R2->glog.AAcnt+R2->glog.TTcnt+R2->glog.GGcnt+R2->glog.CCcnt));
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Jump reads</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mapped-diffChr</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.map_difchr, ((float)OA->glog.map_difchr/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Insertsize>=1K</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.Isize_ge1K, ((float)OA->glog.Isize_ge1K/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Pair orientation</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mapped-pair forward</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.FF, ((float)OA->glog.FF/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mapped-pair reverse</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.RR, ((float)OA->glog.RR/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mapped-pair outward</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.OO, ((float)OA->glog.OO/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mapped-pair innerward</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", (OA->glog.FR-OA->glog.OO), ((float)(OA->glog.FR-OA->glog.OO)/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=mapsummary>\n");
	fprintf(fp,"<h3>Other summary</h3>\n");
	fprintf(fp,"<table class=\"mapsummary infotable\">\n");
	fprintf(fp,"\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mean depth of coverages</td>\n");
	if (prmtr->trgtsize == 0) {fprintf(fp,"<td class=column2>%.2fX</td>\n", (float)(OA->glog.AAcnt+OA->glog.TTcnt+OA->glog.GGcnt+OA->glog.CCcnt)/(float)prmtr->ref_size); }
        else {fprintf(fp,"<td class=column2>%.2fX</td>\n", (float)(OA->glog.AAcnt+OA->glog.TTcnt+OA->glog.GGcnt+OA->glog.CCcnt)/(float)prmtr->trgtsize);}
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>DuplicateMarked rate</td>\n");
	fprintf(fp,"<td class=column2>%ld (%.2f%%)</td>\n", OA->glog.DupCnt, ((float)OA->glog.DupCnt/(float)OA->glog.tot_rd_cnt)*100);
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mean insertsize</td>\n"); 
	fprintf(fp,"<td class=column2>%.1f</td>\n", misize);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mean base quality</td>\n");
	fprintf(fp,"<td class=column2>%.1f</td>\n", mbasq);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mean mapping quality</td>\n");
	fprintf(fp,"<td class=column2>%.1f</td>\n", mmapq);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"<tr onmouseover=\"this.style.backgroundColor='#D5F5E3 ';\" onmouseout=\"this.style.backgroundColor='#FBFCFC';\">\n");
	fprintf(fp,"<td class=column1>Mean read length</td>\n");
	fprintf(fp,"<td class=column2>%.1f</td>\n", mrdln);
	fprintf(fp,"</tr>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</table>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<h2 style=\"background-color:#E5E7E9;\"><u>Summary plots</u></h2>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=\"row\">\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<h3 id=\"Base quality distribution\" style=\"background-color:#FEF9E7;\" onmouseover=\"this.style.backgroundColor='#F9E79F';\" onmouseout=\"this.style.backgroundColor='#FEF9E7';\"> Base quality distribution<a class=\"headerlink\" name=\"Overall_Mismatchplot_U1.png\" title=\"Permalink to this headline\">&nbsp;</a></h3>\n");
	fprintf(fp,"<div><img width=\"500\" height=\"500\" src=\"plots/Basequality.png\"></div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<h3 id=\"Per base average quality\" style=\"background-color:#FEF9E7;\" onmouseover=\"this.style.backgroundColor='#F9E79F';\" onmouseout=\"this.style.backgroundColor='#FEF9E7';\">Per base average quality<a class=\"headerlink\" name=\"Overall_Mismatchplot_U1.png\" title=\"Permalink to this headline\">&nbsp;</a></h3>\n");
	fprintf(fp,"<div><img width=\"500\" height=\"500\" src=\"plots/Percycle_qv.png\"></div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=\"row\">\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<h3 id=\"Per base sequence content\" style=\"background-color:#FEF9E7;\" onmouseover=\"this.style.backgroundColor='#F9E79F';\" onmouseout=\"this.style.backgroundColor='#FEF9E7';\">Per base sequence content<a class=\"headerlink\" name=\"Overall_Mismatchplot_U1.png\" title=\"Permalink to this headline\">&nbsp;</a></h3>\n");
	fprintf(fp,"<div><img width=\"500\" height=\"500\" src=\"plots/Perbase-seq-content.png\"></div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<h3 id=\"GC content distribution\" style=\"background-color:#FEF9E7;\" onmouseover=\"this.style.backgroundColor='#F9E79F';\" onmouseout=\"this.style.backgroundColor='#FEF9E7';\">GC content distribution<a class=\"headerlink\" name=\"Overall_Mismatchplot_U1.png\" title=\"Permalink to this headline\">&nbsp;</a></h3>\n");
	fprintf(fp,"<div><img width=\"500\" height=\"500\" src=\"plots/GC.png\"></div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=\"row\">\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<h3 id=\"Mapping quality distribution\" style=\"background-color:#FEF9E7;\" onmouseover=\"this.style.backgroundColor='#F9E79F';\" onmouseout=\"this.style.backgroundColor='#FEF9E7';\">Mapping quality profile<a class=\"headerlink\" name=\"Overall_Mismatchplot_U1.png\" title=\"Permalink to this headline\">&nbsp;</a></h3>\n");
	fprintf(fp,"<div><img width=\"500\" height=\"500\" src=\"plots/Mapping_quality.png\"></div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<h3 id=\"Insert size histogram\" style=\"background-color:#FEF9E7;\" onmouseover=\"this.style.backgroundColor='#F9E79F';\" onmouseout=\"this.style.backgroundColor='#FEF9E7';\">Insert size histogram<a class=\"headerlink\" name=\"Overall_Mismatchplot_U1.png\" title=\"Permalink to this headline\">&nbsp;</a></h3>\n");
	fprintf(fp,"<div><img width=\"500\" height=\"500\" src=\"plots/Insertsize.png\"></div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<h3 id=\"Mismatch counts\" style=\"background-color:#FEF9E7;\" onmouseover=\"this.style.backgroundColor='#F9E79F';\" onmouseout=\"this.style.backgroundColor='#FEF9E7';\">Mismatch counts<a class=\"headerlink\" name=\"Overall_Mismatchplot_U1.png\" title=\"Permalink to this headline\">&nbsp;</a></h3>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=\"row\">\n");
	fprintf(fp,"  <div class=\"column\">\n");
	fprintf(fp,"      <div><img width=\"600\" height=\"600\" src=\"plots/Overall_Mismatchplot_U.png\"></div>\n");
	fprintf(fp,"  </div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"  <div class=\"column\">\n");
	fprintf(fp,"      <div class=\"row\">\n");
	fprintf(fp,"          <img width=\"300\" height=\"300\" src=\"plots/Read1_fs_Mismatchplot_U.png\">\n");
	fprintf(fp,"          <img width=\"300\" height=\"300\" src=\"plots/Read1_rs_Mismatchplot_U.png\">\n");
	fprintf(fp,"      </div> \n");
	fprintf(fp,"\n");
	fprintf(fp,"      <div class=\"row\">\n");
	fprintf(fp,"          <img width=\"300\" height=\"300\" src=\"plots/Read2_fs_Mismatchplot_U.png\">\n");
	fprintf(fp,"          <img width=\"300\" height=\"300\" src=\"plots/Read2_rs_Mismatchplot_U.png\">\n");
	fprintf(fp,"      </div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"  </div><!-- graph section -->\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<h3 id=\"Basechange & quality\" style=\"background-color:#FEF9E7;\" onmouseover=\"this.style.backgroundColor='#F9E79F';\" onmouseout=\"this.style.backgroundColor='#FEF9E7';\"> Basechange & quality<a class=\"headerlink\" name=\"Overall_Basechange_and_Quality_U1.png\" title=\"Permalink to this headline\">&nbsp;</a></h3>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=\"row\">\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<div><img width=\"600\" height=\"600\" src=\"plots/Overall_Basechange_and_Quality_U.png\"></div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<div class=\"row\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"300\" src=\"plots/Read1_fs_Basechange_and_Quality_U.png\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"300\" src=\"plots/Read1_rs_Basechange_and_Quality_U.png\">\n");
	fprintf(fp,"</div> \n");
	fprintf(fp,"\n");
	fprintf(fp," <div class=\"row\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"300\" src=\"plots/Read2_fs_Basechange_and_Quality_U.png\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"300\" src=\"plots/Read2_rs_Basechange_and_Quality_U.png\">\n");
	fprintf(fp,"  </div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</div><!-- graph section -->  \n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,"<h3 id=\"Per cycle mismatchs & quality status\" style=\"background-color:#FEF9E7;\" onmouseover=\"this.style.backgroundColor='#F9E79F';\" onmouseout=\"this.style.backgroundColor='#FEF9E7';\" > Per cycle mismatchs & quality status<a class=\"headerlink\" name=\"Overall_Com_plots_U.png\" title=\"Permalink to this headline\">&nbsp;</a></h3>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=\"row\">\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<div><img width=\"600\" height=\"550\" src=\"plots/Overall_Com_plots_U.png\"></div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<div class=\"row\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"275\" src=\"plots/Read1_fs_Com_plots_U.png\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"275\" src=\"plots/Read1_rs_Com_plots_U.png\">\n");
	fprintf(fp,"</div> \n");
	fprintf(fp,"\n");
	fprintf(fp," <div class=\"row\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"275\" src=\"plots/Read2_fs_Com_plots_U.png\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"275\" src=\"plots/Read2_rs_Com_plots_U.png\">\n");
	fprintf(fp,"  </div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</div><!-- graph section -->\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<h3 id=\"Deletion and Insertion status\" style=\"background-color:#FEF9E7;\" onmouseover=\"this.style.backgroundColor='#F9E79F';\" onmouseout=\"this.style.backgroundColor='#FEF9E7';\" > Deletion and Insertion status<a class=\"headerlink\" name=\"Overall_Indel_plot_U1.png\" title=\"Permalink to this headline\">&nbsp;</a></h3>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=\"row\">\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<div><img width=\"600\" height=\"600\" src=\"plots/Overall_Indel_plot_U.png\"></div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<div class=\"row\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"300\" src=\"plots/Read1_fs_Indel_plot_U.png\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"300\" src=\"plots/Read1_rs_Indel_plot_U.png\">\n");
	fprintf(fp,"</div> \n");
	fprintf(fp,"\n");
	fprintf(fp," <div class=\"row\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"300\" src=\"plots/Read2_fs_Indel_plot_U.png\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"300\" src=\"plots/Read2_rs_Indel_plot_U.png\">\n");
	fprintf(fp,"  </div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</div><!-- graph section -->\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<h3 id=\"Read clipping status\" style=\"background-color:#FEF9E7;\" onmouseover=\"this.style.backgroundColor='#F9E79F';\" onmouseout=\"this.style.backgroundColor='#FEF9E7';\" > Read clipping status<a class=\"headerlink\" name=\"Overall_Clipped_plot_U1.png\" title=\"Permalink to this headline\">&nbsp;</a></h3>\n");
	fprintf(fp,"\n");
	fprintf(fp,"<div class=\"row\">\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<div><img width=\"600\" height=\"500\" src=\"plots/Overall_Clipped_plot_U.png\"></div>\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"<div class=\"column\">\n");
	fprintf(fp,"<div class=\"row\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"250\" src=\"plots/Read1_fs_Clipped_plot_U.png\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"250\" src=\"plots/Read1_rs_Clipped_plot_U.png\">\n");
	fprintf(fp,"</div> \n");
	fprintf(fp,"\n");
	fprintf(fp," <div class=\"row\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"250\" src=\"plots/Read2_fs_Clipped_plot_U.png\">\n");
	fprintf(fp,"    <img width=\"300\" height=\"250\" src=\"plots/Read2_rs_Clipped_plot_U.png\">\n");
	fprintf(fp,"  </div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"</div><!-- graph section -->\n");
	fprintf(fp,"</div>\n");
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,"    <div class=\"fixed-footer\">\n");
	fprintf(fp,"        <div class=\"container\">NIBMG</div>        \n");
	fprintf(fp,"    </div>\n");
	fprintf(fp,"</body>\n");
	fprintf(fp,"</html>\n");
  fclose(fp);
  free(outfn);
}

 void Arrange(char *outdir)
 {
   char *outfldr=0;
   outfldr = calloc(strlen(outdir) + 500, 1);
   sprintf(outfldr, "mkdir -p %s/plots", outdir);
   system(outfldr);
   *outfldr=0;  
   sprintf(outfldr, "mv %s/*.png %s/plots", outdir, outdir);
   system(outfldr);

   *outfldr=0;
   sprintf(outfldr, "mkdir -p %s/logs", outdir);
   system(outfldr);
   *outfldr=0;  
   sprintf(outfldr, "mv %s/*.txt %s/logs", outdir, outdir);
   system(outfldr);
   free(outfldr);
 }
 void bamqc_help()
 {
   fprintf(stderr, "\n");
   fprintf(stderr, "Program: mapinsights\n");
   fprintf(stderr, "Version: 1.0\n\n");
   fprintf(stderr, "Usage: mapinsights bamqc -r <ref.fa> -o <output-folder-path> -i <aligned.bam>\n\n");
   fprintf(stderr, "Options:\n");
   fprintf(stderr, "	-r  ref.fa	reference fasta\n");        
   fprintf(stderr, "	-b  bed file	regions (BED) [null]\n");
   fprintf(stderr, "	-i  input file	Alignment file (BAM)\n");
   fprintf(stderr, "	-o  path to output folder	[./]\n");
   //fprintf(stderr, "	-t  no of threads	[1]\n");
   fprintf(stderr, "	-x  exclude read groups listed in FILE, one per line [null]\n");
   fprintf(stderr, "	-h  help\n");
   fprintf(stderr, "\n");
   //return 1;
 }

 long int getsiz(char *bname) 
 { 
    FILE* fp = fopen(bname, "r"); 
    
    if (fp == NULL) { 
        printf("File Not Found!\n"); 
        return -1; 
    } 
    fseek(fp, 0L, SEEK_END); 
    long int res = ftell(fp); 
    fclose(fp); 
    return res; 
 } 

 int checkRG(char rg[500][1024], int n, char *id)
 {
    int i=0, f=0;
    //id++;
    while(i < n)
    {
       if(!strcmp(rg[i],id)) {f=1; break;}
       i++;
    }
    return f;
 }
 int main_bamqc(int argc,char *argv[])
 {
   FILE *xfp_rg;
   bamFile fp;
   bam_header_t *xheader;
   bam1_t *b;
   faidx_t *fai;
   read_info_t *ri;
   int tid0 = -1,ref_len,ref_tid = -1, n;
   int tid, ii, Counttt=0;
   int beg, end, progrs_bincnt=1;
   long long int seqcount = 0;
   long int progrs=0;
   //uint32_t *cigar;
   uint64_t mm_cnt1 = 0;
   int x=0, k, rf_pos, rd_pos=0, offset, read_end, xrg_cnt=0;
   char *ref, comnd[500];
   char rg_ch;
   char *ss, *rr, *p, *q;
   int rLen=100, rLenMax=0, read_ovrlp_st, read_ovrlp_nd,i;
   log_v *T,*OA, *RD1, *RD2, *RD1_FS, *RD1_RS, *RD2_FS, *RD2_RS;
   paramtr *prmtr;
   //####### Added on 22Sep2020 ###############
   //log_genrl *glog;
   //long int Acnt=0, Tcnt=0, Gcnt=0, Ccnt=0, Ncnt=0, iSize;
   //long int A1cnt=0, T1cnt=0, G1cnt=0, C1cnt=0, N1cnt=0;
   //char *outdir, *ref_file, *bed_file, *bam_file;
   void *bed = 0;
   
   prmtr = calloc(1, sizeof(paramtr));
   //printf("I am here0\n");
   //prmtr->thrd = 1;
   while ((n = getopt(argc, argv, "r:b:i:o:t:x:h")) >= 0) {
		switch (n) {
			 case 'r': prmtr->ref_file = optarg; break; // reference file
                         case 'b': prmtr->bed_file = optarg; break; // bed file
			 case 'i': prmtr->bam_file = optarg; break; // Input bam file
			 case 'o': prmtr->outdir = optarg; break; // output folder path
                         //case 't': prmtr->thrd = atoi(optarg); break;
                         case 'x': prmtr->xrg_file = optarg; break;
                         case 'h': bamqc_help(); return 1; break;
                        }
	}
	if (argc == 1 || optind-1 == argc) {bamqc_help(); return 1;}
   
   if(access(prmtr->ref_file, F_OK) == -1) {printf("Reference file does not exist.\n"); return 1;}
      
   if(access(prmtr->bam_file, F_OK) == -1) {printf("Alignment(.bam) file does not exist.\n"); return 1;}
   
   if(prmtr->bed_file) {if(access(prmtr->bed_file, F_OK) == -1) {printf("Bed does not exist.\n"); return 1;}}

   if(access(prmtr->outdir, F_OK) == -1) {printf("Output filder does not exist.\n"); return 1;}
   
   if(prmtr->xrg_file) {if(access(prmtr->xrg_file, F_OK) == -1) {printf("Read-groups file does not exist.\n"); return 1;}}
   
   //if(prmtr->thrd)
   ii=0;
   for(i=0;i<argc; i++) { /*printf("%s ",argv[i]);*/ ii+=strlen(argv[i]); }
   printf("\n");
   *comnd=0;
   prmtr->comnd = calloc(ii + 50, 1);
   for(i=0;i<argc; i++) { strcat(prmtr->comnd, argv[i]); strcat(prmtr->comnd, " ");}

   printf("Program: mapinsights\n");
   printf("Version: 1.0\n\n");
   printf("command :: mapinsights %s\n\n", prmtr->comnd);

   //if(prmtr->thrd < 1) {printf("%d thread not allowed.\n", prmtr->thrd); return 1;}
   //xcore=bam_init1();
   //printf("I am here\n");
   //printf("No of threads = %d\n", prmtr->thrd);
   //ref_size = 3137454505;
   seqcount = getseqcount(prmtr->bam_file);
   
   fp = bam_open(prmtr->bam_file, "r");
   //printf("I am here1\n");
   if(fp==NULL) { printf("Unable to open BAM file\n"); return 1;}
   fai = fai_load(prmtr->ref_file);
   if (fai == 0) { printf("Unable to open reference file\n"); return 1;}
   //printf("I am here2\n");
   if(prmtr->bed_file) {bed = bed_read_new(prmtr->bed_file, &prmtr->trgtsize);} //BED
   
   ii=0;
   if(prmtr->xrg_file) 
   {
     xfp_rg=fopen(prmtr->xrg_file, "r");
     if(xfp_rg == NULL) {printf("Unable to open read-group file\n"); return 1;}
     prmtr->xrg_cnt=0; k=0;
     while((rg_ch=fgetc(xfp_rg))!=EOF)
     {
       if(rg_ch=='\n') {prmtr->x_rg[xrg_cnt][k]='\0'; k=0; prmtr->xrg_cnt++;}
       prmtr->x_rg[prmtr->xrg_cnt][k]=rg_ch; k++;
     }
     prmtr->x_rg[prmtr->xrg_cnt][k]='\0'; prmtr->xrg_cnt++;
     printf("RG exclude list\n");
     for(k=0; k<prmtr->xrg_cnt; k++) {printf("%s\n",prmtr->x_rg[k]);}
   }
   //printf("I am here3\n");
      beg = 0; end = 1<<30; tid = -1;
      xheader=bam_header_init();
      xheader=bam_header_read(fp);
      for(i=0; i< xheader->n_targets; i++)
      {
        prmtr->ref_size+=xheader->target_len[i];
        //printf("%s	%d\n", xheader->target_name[i], xheader->target_len[i]);
      }
      prmtr->no_of_contigs = xheader->n_targets;
      b=bam_init1();
      //printf("%s\n",xheader->text);
      if ((rr = strstr(xheader->text, "@RG\t")) == 0) 
      {
          printf("RG tag is missing\n");
      }
      else 
      {
      
        if ((ss = strstr(rr, "\tSM:")) != 0) {ss += 4;
        k=0;
        while(ss[k]) {if(ss[k] == '\t' || ss[k] == '\n') {break;} k++;}
        p=calloc(k+100,1);
        strncpy(p, ss, k);
        prmtr->sample_name=p;}
        
        
        if ((ss = strstr(rr, "\tLB:")) != 0) {ss += 4;
        k=0;
        while(ss[k]) {if(ss[k] == '\t' || ss[k] == '\n') {break;} k++;}
        q=calloc(k+50,1);
        strncpy(q, ss, k);
        prmtr->library=q;}
        
        while((rr = strstr(rr, "@RG\t")) !=0) 
        {
           if((ss = strstr(rr, "ID:")) !=0)
           {
                ss+=3;k=0;
                while(*ss!='\t' && *ss!='\n')
                {
                   prmtr->sId[prmtr->no_of_rgs][k]=*ss;
                   ss++;k++;
                }
                prmtr->sId[prmtr->no_of_rgs][k]='\0';
           }
           rr=ss+k; prmtr->no_of_rgs++;
        }
        if(k>=1000) {printf("@RG ID: tag is too big\n"); return 1;} else {prmtr->sId[prmtr->no_of_rgs][k]='\0';}
        if(prmtr->no_of_rgs>=499) {printf("Too many @RG tags\n"); return 1;}
         
        
        
      }
      
      //printf("Library =%s\n",prmtr->library);
      //printf("Sample nmae =%s\n",prmtr->sample_name);
      //printf("No of rg = %d\n", prmtr->no_of_rgs);
      //for(k=0;k<prmtr->no_of_rgs;k++)
      //{printf("%s\n",prmtr->sId[k]);}
      //*ss='\0';
      //printf("sample length = %d\n",k);
      
      ri = calloc(1, sizeof(read_info_t));
      T = calloc(1, sizeof(log_v));
      OA = calloc(1, sizeof(log_v));
      RD1 = calloc(1, sizeof(log_v));
      RD2 = calloc(1, sizeof(log_v));
      RD1_FS = calloc(1, sizeof(log_v));
      RD1_RS = calloc(1, sizeof(log_v));
      RD2_FS = calloc(1, sizeof(log_v));
      RD2_RS = calloc(1, sizeof(log_v));
      
         //printf("Target size = %ld\n", trgtsize);
        tid0 = tid;
        if (tid0 >= 0 && fai) { // region is set
		ref = faidx_fetch_seq(fai, xheader->target_name[tid0], 0, 0x7fffffff, &ref_len);
		ref_tid = tid0;
		
	} else ref_tid = -1, ref = 0;
      printf("QC analysis is in progress....\n");
      while( (x=bam_read1(fp,b)) >= 0 )
      {
         
         uint8_t *seq = bam1_seq(b); //, *qual = bam1_qual(b);
         
         uint32_t *cigar = bam_get_cigar(b);
         //cumsize += x; cumsize1+=x-4;
         //printf("x= %d	x-4= %d cumsize= %ld\n", x, x-4, cumsize);
         //printf("%.2f %% complete\n", ((float)OA->glog.all_rd_cnt/(float)seqcount)*100);
         
         if( progrs == (int)(5*seqcount/100)+5) {printf("%d%% complete\n",5*progrs_bincnt); progrs=0; progrs_bincnt++;}
         progrs++;
         OA->glog.all_rd_cnt++; 
         if (b->core.flag & BAM_FUNMAP) {OA->glog.unmap++; continue;}
         
         read_ovrlp_st = 0; read_ovrlp_nd = b->core.l_qseq;
         //printf("Initial::%d     %d\n",read_ovrlp_st, read_ovrlp_nd);
         read_end = bam_calend(&b->core, bam1_cigar(b));
         if (bed && bed_overlap_u(bed, xheader->target_name[b->core.tid], b->core.pos, read_end, &read_ovrlp_st, &read_ovrlp_nd) == 0) continue;
         //printf("Read group = %s\n",bam_aux_get(b, "RG"));
         //printf("rdst=%d rdnd=%d\n", read_ovrlp_st, read_ovrlp_nd);
         if(prmtr->xrg_file && prmtr->xrg_cnt > 0) {
         ii=0; 
         uint8_t *rg=bam_aux_get(b, "RG");
         if(rg) { 
         //while(rr[ii+1] != '\0') { r_rg[ii]=rr[ii+1]; ii++;} r_rg[ii]='\0';
         rg++;
         if(checkRG(prmtr->x_rg, prmtr->xrg_cnt, rg)) { printf("Exclude SAM records with RG:ID = %s\n", rg); continue;}
         }
         }
         
         offset = (((read_end - b->core.pos)) - b->core.l_qseq);
         
         if(read_ovrlp_st < 0 ) 
         {
            read_ovrlp_st = read_ovrlp_st + offset; 
            if(read_ovrlp_st < 0 ) {read_ovrlp_st = 0;}
         }
         else if(read_ovrlp_st > b->core.l_qseq) 
         {
           read_ovrlp_st = (read_ovrlp_st - offset);
           if(read_ovrlp_st > b->core.l_qseq) {read_ovrlp_st = b->core.l_qseq;}
         }
         else if(read_ovrlp_nd < 0 ) 
         {
           read_ovrlp_nd = read_ovrlp_nd + offset;
           if(read_ovrlp_nd < 0 ) {read_ovrlp_nd = 0;} 
         }
         else if(read_ovrlp_nd > b->core.l_qseq) 
         {
           read_ovrlp_nd = read_ovrlp_nd - offset;
           if(read_ovrlp_nd > b->core.l_qseq) {read_ovrlp_nd = b->core.l_qseq;}
         }
         

         Counttt++;
	 rLen=b->core.l_qseq;
         if(rLenMax < rLen) {rLenMax = rLen;}
         
         if (b->core.tid != ref_tid) {
			free(ref); ref = 0;
                       if (fai) ref = faidx_fetch_seq(fai, xheader->target_name[b->core.tid], 0, 0x7fffffff, &ref_len);
			ref_tid = b->core.tid;
                        
		}
         
         ri->dif_count=0;
         //////////Test mismatch counts/////////
         rd_pos=0;rf_pos=b->core.pos;
          //printf("I am here\n");
          for (k = 0 ; k < b->core.n_cigar ; ++k)
          {

                int l = cigar[k]>>4;
                int op = cigar[k]&0xf;
                //printf("I am here1\n");
                if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
                {
                   
                   for( i=0; i<l; i++)
                   {
                       
                       if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)] != ((ref && (rf_pos+i) < ref_len)? toupper(ref[rf_pos+i]) : 'N'))
                       {
                       
                         if(bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)]!='N' || toupper(ref[rf_pos+i])!='N')
                         {
                          //printf("%s	%d	%c %c %c -> %c %c %c\n",xheader->target_name[b->core.tid], rf_pos+i, toupper(ref[rf_pos+i-1]), toupper(ref[rf_pos+i]), toupper(ref[rf_pos+i+1]), toupper(ref[rf_pos+i-1]), bam_nt16_rev_table[bam1_seqi(seq,rd_pos+i)], toupper(ref[rf_pos+i+1]));
                          mm_cnt1+=1;}
                       }
                   }
                    rf_pos+=l; rd_pos+=l;
                }
                else if(op == BAM_CDEL || op == BAM_CREF_SKIP) {rf_pos+=l;}
                else if(op == BAM_CINS) {rd_pos+=l;}
                else if(op == BAM_CSOFT_CLIP){rd_pos+=l;}
                 //printf("I am here2\n");

         }

          //printf("I am here3\n");
         ///////////////////////////////////////
         if(ref_tid==b->core.tid && ref && fai){
         OA->cnt+=1;
         get_seqInfo(b, ref, ref_len, T, OA, read_ovrlp_st, read_ovrlp_nd);
         //Addup_u(OA, T, 0, rLen, "Overall");
         //printf("Mapping quality=%d	Zero=%ld	1-10=%ld	Sixty=%ld\n",b->core.qual, OA->glog.MapQ_bin_dist[1], OA->glog.MapQ_bin_dist[2], OA->glog.MapQ_bin_dist[7]);
         if(b->core.flag & BAM_FREAD1) 
         {
            //Addup_u(RD1, T, 0, rLen, "Read1");
            if(b->core.flag & BAM_FREVERSE) 
            {
                Addup_u_u(OA, RD1, RD1_RS, T, 1, rLen, 1);
                RD1_RS->cnt+=1;
                //Addup_u(RD1_RS, T, 1, rLen, "Read1_rs");
            }
            else
            {
                Addup_u_u(OA, RD1, RD1_FS, T, 1, rLen, 0);
                RD1_FS->cnt+=1;
                //Addup_u(RD1_FS, T, 1, rLen, "Read1_fs");            
            }

         }
         else
         {
            //Addup_u(RD2, T, 0, rLen, "Read2");
            if(b->core.flag & BAM_FREVERSE) 
            {
                Addup_u_u(OA, RD2, RD2_RS, T, 1, rLen, 1);
                RD2_RS->cnt+=1;
                //Addup_u(RD2_RS, T, 1, rLen, "Read2_rs");
            }
            else
            {
                Addup_u_u(OA, RD2, RD2_FS, T, 1, rLen, 0);
                RD2_FS->cnt+=1;
                //Addup_u(RD2_FS, T, 1, rLen, "Read2_fs");            
            } 
         }
         
         /*get_seqInfo(b, ref, ref_len , ri);*/}
         
         //printf("Miss-match = %d	%d	%d	%d	%s	",ri->dif_count, ri->sq_len, ri->sc_pct, ri->q_sum,bam_aux_get(b, "MD"));
         //ch=bam_format1(xheader,b);
         //puts(ch);
         
      }
      bam_close(fp);
      
      //DisLog(OA);
      if(bed) {PrintGenrlLogFiles_u(OA, rLenMax, prmtr->trgtsize, prmtr->outdir, "Overall", prmtr);} else {PrintGenrlLogFiles_u(OA, rLenMax, prmtr->ref_size, prmtr->outdir, "Overall", prmtr);}
      PrintGenrlLogFiles_u(RD1, rLenMax, prmtr->trgtsize, prmtr->outdir, "Read1", prmtr);
      PrintGenrlLogFiles_u(RD2, rLenMax, prmtr->trgtsize, prmtr->outdir, "Read2", prmtr);
      PrintLogFiles(OA, 1, rLenMax, prmtr->outdir);
      makeplots_v1(prmtr->outdir, "Overall");
      //printf("Mean depth of coverage = %f\n",(float)(T->glog.AAcnt+T->glog.TTcnt+T->glog.GGcnt+T->glog.CCcnt)/(float)trgtsize);
      PrintLogFiles(RD1_FS, 2, rLenMax, prmtr->outdir);
      makeplots_v1(prmtr->outdir, "Read1_fs");
      PrintLogFiles(RD1_RS, 3, rLenMax, prmtr->outdir);
      makeplots_v1(prmtr->outdir, "Read1_rs");
      PrintLogFiles(RD2_FS, 4, rLenMax, prmtr->outdir);
      makeplots_v1(prmtr->outdir, "Read2_fs");
      PrintLogFiles(RD2_RS, 5, rLenMax, prmtr->outdir);
      makeplots_v1(prmtr->outdir, "Read2_rs");
      //system("R CMD BATCH Mismatch_plot_v2_All.R");
      //printf("Total no of reads =%d\n",Counttt);
      //printf("No of mismatch = %ld\n",mm_cnt1);
      //printf("A counts = %ld\n", Acnt);
      //printf("T counts = %ld\n", Tcnt);
      //printf("G counts = %ld\n", Gcnt);
      //printf("C counts = %ld\n", Ccnt);
      //printf("N counts = %ld\n", Ncnt);

      //printf("A1 counts = %ld\n", A1cnt);
      //printf("T1 counts = %ld\n", T1cnt);
      //printf("G1 counts = %ld\n", G1cnt);
      //printf("C1 counts = %ld\n", C1cnt);
      //printf("N1 counts = %ld\n", N1cnt);
  
      

  
      
      
      
      //printf("Level1\n");
      draw_plots(prmtr->outdir);
      //printf("Level2\n");
      Arrange(prmtr->outdir);
      //printf("Level3\n");
      
      time_t t;
      time(&t);
      //printf("Level3.1\n");
      crthtml(prmtr, ctime(&t), OA, RD1, RD2);
      //printf("Analysis date & time = %s\n", ctime(&t));
      
      //printf("Total no of mismatchs = %d\n", T->midc.basemmc);
      //for(x = 0 ; x < 12 ; x++) {printf("%d\n",T->MMC[x]);}
      
      if (bed) bed_destroy_u(bed);
      //if(bed) {printf("Bed is allocated\n");}
      //printf("Level4\n");
      //printf("SeqCounts = %llu\n", seqcount);
      //printf("line Count = %ld\n", OA->glog.all_rd_cnt);
      printf("100%% complete\n\n");
      printf("QC analysis done\n\n");
      return 0;
           
 }

