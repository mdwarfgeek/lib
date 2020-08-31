# Makefile for liblfa

### Constants: edit to suit your system ###

# C compiler
#CC=gcc

# Language flags
LANG=-std=gnu99 -fno-strict-aliasing

# Optimization flags

# DEBUG:
#OPT=-g

# OPT:
OPT=-g -O3

# Compiler flags
CFLAGS=$(LANG) -Wall $(OPT) -fPIC -I$(HOME)/include

# Extra flags for CFITSIO
CFITSIO_INC?=`pkg-config cfitsio --cflags`
CFITSIO_LIBS?=`pkg-config cfitsio --libs`

# FFTW3 float version include and library paths
FFTW3F_INC?=`pkg-config fftw3f --cflags`
FFTW3F_LIBS?=`pkg-config fftw3f --libs`

# Termcap libraries.  Set blank on Win32.
TERMCAP_LIBS?=-lncurses

### End constants section ###

SRCS=airmass.c ap.c background.c bary.c cholesky.c cnf.c const.c cvtunit.c dict.c dplate.c dtai.c dtdb.c filt1d.c fpcoord.c fsgp.c geoc.c iers.c jpleph.c kepler.c matvec.c mjd.c mount.c mpc.c observer.c parallactic.c pixovcirc.c prenut.c qr.c refract.c rng.c shuffle.c skylevel.c sort.c source.c specfind.c strutil.c stumpff.c sysinfo.c util.c wcs.c
OBJS=${SRCS:%.c=%.o}

EXTRA_SRCS=fftutil.c fitsutil.c lautil.c tcutil.c
EXTRA_OBJS=${EXTRA_SRCS:%.c=%.o}

TESTAP_SRCS=testap.c fitsutil.c
TESTAP_OBJS=${TESTAP_SRCS:%.c=%.o}

TESTBACKGROUND_SRCS=testbackground.c fitsutil.c
TESTBACKGROUND_OBJS=${TESTBACKGROUND_SRCS:%.c=%.o}

TESTMPC_SRCS=testmpc.c tcutil.c
TESTMPC_OBJS=${TESTMPC_SRCS:%.c=%.o}

TESTPLAN_SRCS=testplan.c tcutil.c
TESTPLAN_OBJS=${TESTPLAN_SRCS:%.c=%.o}

TESTRNG_SRCS=testrng.c
TESTRNG_OBJS=${TESTRNG_SRCS:%.c=%.o}

TESTSIMPLE_SRCS=testsimple.c
TESTSIMPLE_OBJS=${TESTSIMPLE_SRCS:%.c=%.o}

# Rules for building

.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

all: liblfa.a

depend:
	$(CC) $(CFLAGS) -E -MM $(SRCS) > .depend

extra: $(EXTRA_OBJS)

liblfa.a: $(OBJS)
	rm -f $@
	ar r $@ $(OBJS)
	ranlib $@

fftutil.o: fftutil.c
	$(CC) $(CFLAGS) $(FFTW3F_INC) -o $@ -c $<

fitsutil.o: fitsutil.c
	$(CC) $(CFLAGS) $(CFITSIO_INC) -o $@ -c $<

testap: $(TESTAP_OBJS) liblfa.a
	$(CC) -o testap $(TESTAP_OBJS) liblfa.a $(CFITSIO_LIBS) -lm

testbackground: $(TESTBACKGROUND_OBJS) liblfa.a
	$(CC) -o testbackground $(TESTBACKGROUND_OBJS) \
                 liblfa.a $(CFITSIO_LIBS) -lm

testmpc: $(TESTMPC_OBJS) liblfa.a
	$(CC) -o testmpc $(TESTMPC_OBJS) liblfa.a $(TERMCAP_LIBS) -lm

testplan: $(TESTPLAN_OBJS) liblfa.a
	$(CC) -o testplan $(TESTPLAN_OBJS) liblfa.a $(TERMCAP_LIBS) -lm

testrng: $(TESTRNG_OBJS) liblfa.a
	$(CC) -o testrng $(TESTRNG_OBJS) liblfa.a -lm

testsimple: $(TESTSIMPLE_OBJS) liblfa.a
	$(CC) -o testsimple $(TESTSIMPLE_OBJS) liblfa.a -lm

clean:
	rm -f $(OBJS) liblfa.a
	rm -f $(TESTAP_OBJS) testap
	rm -f $(TESTBACKGROUND_OBJS) testbackground
	rm -f $(TESTMPC_OBJS) testmpc
	rm -f $(TESTPLAN_OBJS) testplan
	rm -f $(TESTRNG_OBJS) testrng
	rm -f $(TESTSIMPLE_OBJS) testsimple
	rm -f $(EXTRA_OBJS)
	rm -f .depend
