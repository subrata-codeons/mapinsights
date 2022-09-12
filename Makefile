CC=			gcc
CFLAGS=		-g -O2 #-Wall
#LDFLAGS=		-Wl,-rpath,\$$ORIGIN/../lib
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_CURSES_LIB=0

LOBJS=		bgzf.o kstring.o bam_aux.o bam.o bam_import.o sam.o bam_index.o	\
			bam_pileup.o razf.o faidx.o \
			sam_header.o
AOBJS=		bamqc.o allmod.o genedepth.o siteinfo.o \
		jumpreads.o multisamplebamqc.o batchplotbamqc.o
PROG=		mapinsights
INCLUDES=	-I.
SUBDIRS=	.
LIBPATH=
#LIBCURSES=	-lcurses # -lXCurses

.SUFFIXES:.c .o
.PHONY: all lib

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all-recur lib-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" LIBPATH="$(LIBPATH)" $$target || exit 1; \
			cd $$wdir; \
		done;

all:$(PROG)

.PHONY:all lib clean cleanlocal
.PHONY:all-recur lib-recur clean-recur cleanlocal-recur install-recur

lib:libbam.a

libbam.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

mapinsights:lib-recur $(AOBJS)
		$(CC) $(CFLAGS) -o $@ $(AOBJS) $(LDFLAGS) libbam.a $(LIBPATH) $(LIBCURSES) -lm -lz -lpthread


bgzf.o:bgzf.c bgzf.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -DBGZF_CACHE $(INCLUDES) bgzf.c -o $@

razf.o:razf.h		
bam.o:bam.h razf.h bam_endian.h kstring.h sam_header.h
sam.o:sam.h bam.h
kstring.o:kstring.h
bam_import.o:bam.h kseq.h khash.h razf.h
bam_pileup.o:bam.h razf.h ksort.h
bam_index.o:bam.h khash.h ksort.h razf.h bam_endian.h
sam_header.o:sam_header.h khash.h
faidx.o:faidx.h razf.h khash.h
bam_aux.o:bam.h khash.h
bamqc.o:bam.h faidx.h kstring.h ksort.h kseq.h khash.h
genedepth.o:bam.h faidx.h
jumpreads.o:bam.h
siteinfo.o:bam.h faidx.h
multisamplebamqc.o:bam.h
batchplotbamqc:bam.h

libbam.1.dylib-local:$(LOBJS)
		libtool -dynamic $(LOBJS) -o libbam.1.dylib -lc -lz

libbam.so.1-local:$(LOBJS)
		$(CC) -shared -Wl,-soname,libbam.so -o libbam.so.1 $(LOBJS) -lc -lz

dylib:
		@$(MAKE) cleanlocal; \
		case `uname` in \
			Linux) $(MAKE) CFLAGS="$(CFLAGS) -fPIC" libbam.so.1-local;; \
			Darwin) $(MAKE) CFLAGS="$(CFLAGS) -fPIC" libbam.1.dylib-local;; \
			*) echo 'Unknown OS';; \
		esac


cleanlocal:
		rm -fr gmon.out *.o a.out *.exe *.dSYM razip bgzip $(PROG) *~ *.a *.so.* *.so *.dylib

clean:cleanlocal-recur


