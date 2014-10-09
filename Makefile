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
CFITSIO_INC?=-I/usr/local/include

# Termcap libraries.  Set blank on Win32.
TERMCAP_LIBS?=-lncurses

### End constants section ###

SRCS=airmass.c bary.c cholesky.c cnf.c cvtunit.c dplate.c dtai.c dtdb.c fpcoord.c geoc.c iers.c jpleph.c linsolve.c matvec.c mjd.c mount.c mpc.c observer.c parallactic.c prenut.c refract.c rng.c shuffle.c skylevel.c sort.c source.c strutil.c stumpff.c sysinfo.c util.c wcs.c
OBJS=${SRCS:%.c=%.o}

EXTRA_SRCS=fitsutil.c lautil.c tcutil.c
EXTRA_OBJS=${EXTRA_SRCS:%.c=%.o}

TESTMPC_SRCS=testmpc.c tcutil.c
TESTMPC_OBJS=${TESTMPC_SRCS:%.c=%.o}

TESTPLAN_SRCS=testplan.c tcutil.c
TESTPLAN_OBJS=${TESTPLAN_SRCS:%.c=%.o}

TESTSIMPLE_SRCS=testsimple.c
TESTSIMPLE_OBJS=${TESTSIMPLE_SRCS:%.c=%.o}

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
	$(CC) $(CFLAGS) $(CFITSIO_INC) -o $@ -c $<

testmpc: $(TESTMPC_OBJS) liblfa.a
	$(CC) -o testmpc $(TESTMPC_OBJS) liblfa.a $(TERMCAP_LIBS) -lm

testplan: $(TESTPLAN_OBJS) liblfa.a
	$(CC) -o testplan $(TESTPLAN_OBJS) liblfa.a $(TERMCAP_LIBS) -lm

testsimple: $(TESTSIMPLE_OBJS) liblfa.a
	$(CC) -o testsimple $(TESTSIMPLE_OBJS) liblfa.a -lm

clean:
	rm -f $(OBJS) liblfa.a
	rm -f $(TESTMPC_OBJS) testmpc
	rm -f $(TESTPLAN_OBJS) testplan
	rm -f $(TESTSIMPLE_OBJS) testsimple
	rm -f $(EXTRA_OBJS)
	rm -f .depend

sinclude .depend
