# Makefile for liblfa: lightweight fundamental astronomy library.

### Constants: edit to suit your system ###

# C compiler
#CC=gcc

# Compiler flags
CFLAGS=-Wall -Wno-strict-aliasing -g -O3 -ffast-math -fPIC -I$(HOME)/include

# Recommended flags for modern cpu
# DEBUG: -g
# OPTS: -g -O3 -ffast-math

# Extra flags for CFITSIO
CFITSIO_CFLAGS=-I/usr/local/include

### End constants section ###

SRCS=airmass.c bary.c dtai.c dtdb.c dsincos.c fpcoord.c geoc.c iers.c jpleph.c matvec.c mjd.c mount.c mpc.c observer.c parallactic.c prenut.c refract.c rng.c skylevel.c sort.c source.c strutil.c stumpff.c wcs.c
OBJS=${SRCS:%.c=%.o}

EXTRA_SRCS=cnf.c cvtunit.c fitsutil.c tcutil.c util.c
EXTRA_OBJS=${EXTRA_SRCS:%.c=%.o}

TESTDSINCOS_SRCS=testdsincos.c
TESTDSINCOS_OBJS=${TESTDSINCOS_SRCS:%.c=%.o}

TESTMPC_SRCS=testmpc.c cvtunit.c util.c
TESTMPC_OBJS=${TESTMPC_SRCS:%.c=%.o}

TESTJPL_SRCS=testjpl.c cvtunit.c util.c
TESTJPL_OBJS=${TESTJPL_SRCS:%.c=%.o}

TESTOBS_SRCS=testobs.c cvtunit.c util.c
TESTOBS_OBJS=${TESTOBS_SRCS:%.c=%.o}

TESTPLAN_SRCS=testplan.c cvtunit.c tcutil.c util.c
TESTPLAN_OBJS=${TESTPLAN_SRCS:%.c=%.o}

TESTREFRO_SRCS=testrefro.c refract.c util.c
TESTREFRO_OBJS=${TESTREFRO_SRCS:%.c=%.o}

TESTSTUMPFF_SRCS=teststumpff.c util.c
TESTSTUMPFF_OBJS=${TESTSTUMPFF_SRCS:%.c=%.o}

TESTSUN_SRCS=testsun.c cvtunit.c util.c
TESTSUN_OBJS=${TESTSUN_SRCS:%.c=%.o}

TESTTP_SRCS=testtp.c cvtunit.c util.c
TESTTP_OBJS=${TESTTP_SRCS:%.c=%.o}

# Rules for building

.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

all: liblfa.a

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
	$(CC) -o testmpc $(TESTMPC_OBJS) -L$(HOME)/lib64 liblfa.a -lsla -lm

testjpl: $(TESTJPL_OBJS) liblfa.a
	$(CC) -o testjpl $(TESTJPL_OBJS) -L$(HOME)/lib64 liblfa.a -lsofa_c -lm

testobs: $(TESTOBS_OBJS) liblfa.a
	$(CC) -o testobs $(TESTOBS_OBJS) -L$(HOME)/lib64 liblfa.a -lsla -lm

testplan: $(TESTPLAN_OBJS) liblfa.a
	$(CC) -o testplan $(TESTPLAN_OBJS) liblfa.a -lncurses -lm

testrefro: $(TESTREFRO_OBJS) liblfa.a
	$(CC) -o testrefro $(TESTREFRO_OBJS) -L$(HOME)/lib64 liblfa.a -lsla -lm

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
	rm -f $(TESTSUN_OBJS) testsun
	rm -f $(TESTSTUMPFF_OBJS) teststumpff
	rm -f $(TESTTP_OBJS) testtp
	rm -f $(EXTRA_OBJS)
