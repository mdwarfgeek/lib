# Makefile for liblfa

### Constants: edit to suit your system ###

# C compiler
#CC=gcc

# Optimization flags

# DEBUG:
#OPT=-g

# OPT:
OPT=-g -O3 -ffast-math

# Compiler flags
CFLAGS=-std=gnu99 -Wall -Wno-strict-aliasing $(OPT) -fPIC -I$(HOME)/include

# Extra flags for CFITSIO
CFITSIO_CFLAGS=-I/usr/local/include

### End constants section ###

SRCS=airmass.c bary.c cnf.c cvtunit.c dtai.c dtdb.c dsincos.c fpcoord.c geoc.c iers.c jpleph.c matvec.c mjd.c mount.c mpc.c observer.c parallactic.c prenut.c refract.c rng.c skylevel.c sort.c source.c strutil.c stumpff.c sysinfo.c util.c wcs.c
OBJS=${SRCS:%.c=%.o}

EXTRA_SRCS=fitsutil.c lautil.c tcutil.c
EXTRA_OBJS=${EXTRA_SRCS:%.c=%.o}

TESTDSINCOS_SRCS=testdsincos.c
TESTDSINCOS_OBJS=${TESTDSINCOS_SRCS:%.c=%.o}

TESTMPC_SRCS=testmpc.c tcutil.c
TESTMPC_OBJS=${TESTMPC_SRCS:%.c=%.o}

TESTJPL_SRCS=testjpl.c
TESTJPL_OBJS=${TESTJPL_SRCS:%.c=%.o}

TESTOBS_SRCS=testobs.c
TESTOBS_OBJS=${TESTOBS_SRCS:%.c=%.o}

TESTPLAN_SRCS=testplan.c tcutil.c
TESTPLAN_OBJS=${TESTPLAN_SRCS:%.c=%.o}

TESTREFRO_SRCS=testrefro.c
TESTREFRO_OBJS=${TESTREFRO_SRCS:%.c=%.o}

TESTSIMPLE_SRCS=testsimple.c
TESTSIMPLE_OBJS=${TESTSIMPLE_SRCS:%.c=%.o}

TESTSTUMPFF_SRCS=teststumpff.c
TESTSTUMPFF_OBJS=${TESTSTUMPFF_SRCS:%.c=%.o}

TESTSUN_SRCS=testsun.c
TESTSUN_OBJS=${TESTSUN_SRCS:%.c=%.o}

TESTTP_SRCS=testtp.c
TESTTP_OBJS=${TESTTP_SRCS:%.c=%.o}

# Rules for building

.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

all: liblfa.a

depend:
	$(CC) -E -MM $(SRCS) > .depend

extra: $(EXTRA_OBJS)

liblfa.a: $(OBJS)
	rm -f $@
	ar r $@ $(OBJS)
	ranlib $@

fitsutil.o: fitsutil.c
	$(CC) $(CFLAGS) $(CFITSIO_CFLAGS) -o $@ -c $<

testdsincos: $(TESTDSINCOS_OBJS) liblfa.a
	$(CC) -o testdsincos $(TESTDSINCOS_OBJS) liblfa.a -lm

testmpc: $(TESTMPC_OBJS) liblfa.a
	$(CC) -o testmpc $(TESTMPC_OBJS) -L$(HOME)/lib64 liblfa.a -lsla -lncurses -lm

testjpl: $(TESTJPL_OBJS) liblfa.a
	$(CC) -o testjpl $(TESTJPL_OBJS) -L$(HOME)/lib64 liblfa.a -lsofa_c -lm

testobs: $(TESTOBS_OBJS) liblfa.a
	$(CC) -o testobs $(TESTOBS_OBJS) -L$(HOME)/lib64 liblfa.a -lsla -lm

testplan: $(TESTPLAN_OBJS) liblfa.a
	$(CC) -o testplan $(TESTPLAN_OBJS) liblfa.a -lncurses -lm

testrefro: $(TESTREFRO_OBJS) liblfa.a
	$(CC) -o testrefro $(TESTREFRO_OBJS) -L$(HOME)/lib64 liblfa.a -lsla -lm

testsimple: $(TESTSIMPLE_OBJS) liblfa.a
	$(CC) -o testsimple $(TESTSIMPLE_OBJS) liblfa.a -lm

teststumpff: $(TESTSTUMPFF_OBJS) liblfa.a
	$(CC) -o teststumpff $(TESTSTUMPFF_OBJS) -L$(HOME)/lib64 liblfa.a -lsla -lm

testsun: $(TESTSUN_OBJS) liblfa.a
	$(CC) -o testsun $(TESTSUN_OBJS) -L$(HOME)/lib64 liblfa.a -lsla -lm

testtp: $(TESTTP_OBJS) liblfa.a
	$(CC) -o testtp $(TESTTP_OBJS) -L$(HOME)/lib64 liblfa.a -lsla -lm

clean:
	rm -f $(OBJS) liblfa.a
	rm -f $(TESTDSINCOS_OBJS) testdsincos
	rm -f $(TESTMPC_OBJS) testmpc
	rm -f $(TESTJPL_OBJS) testjpl
	rm -f $(TESTOBS_OBJS) testobs
	rm -f $(TESTPLAN_OBJS) testplan
	rm -f $(TESTREFRO_OBJS) testrefro
	rm -f $(TESTSIMPLE_OBJS) testsimple
	rm -f $(TESTSUN_OBJS) testsun
	rm -f $(TESTSTUMPFF_OBJS) teststumpff
	rm -f $(TESTTP_OBJS) testtp
	rm -f $(EXTRA_OBJS)
	rm -f .depend

sinclude .depend
