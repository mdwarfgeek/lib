# Makefile for liblfa: lightweight fundamental astronomy library.

### Constants: edit to suit your system ###

# C compiler
#CC=gcc

# Compiler flags
CFLAGS=-Wall -g -fPIC -I$(HOME)/include

# Recommended flags for modern cpu
# DEBUG: -g
# OPTS: -O3 -msse3

### End constants section ###

SRCS=airmass.c bary.c dtai.c dtdb.c geoc.c iers.c jpleph.c matvec.c observer.c prenut.c refract.c source.c strutil.c
OBJS=${SRCS:%.c=%.o}

TESTJPL_SRCS=testjpl.c cvtunit.c util.c
TESTJPL_OBJS=${TESTJPL_SRCS:%.c=%.o}

TESTOBS_SRCS=testobs.c cvtunit.c util.c
TESTOBS_OBJS=${TESTOBS_SRCS:%.c=%.o}

TESTPLAN_SRCS=testplan.c cvtunit.c tcutil.c util.c
TESTPLAN_OBJS=${TESTPLAN_SRCS:%.c=%.o}

TESTREFRO_SRCS=testrefro.c refract.c util.c
TESTREFRO_OBJS=${TESTREFRO_SRCS:%.c=%.o}

TESTSUN_SRCS=testsun.c cvtunit.c util.c
TESTSUN_OBJS=${TESTSUN_SRCS:%.c=%.o}

# Rules for building

.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

all: liblfa.a

liblfa.a: $(OBJS)
	rm -f $@
	ar r $@ $(OBJS)
	ranlib $@

testjpl: $(TESTJPL_OBJS) liblfa.a
	$(CC) -o testjpl $(TESTJPL_OBJS) -L$(HOME)/lib64 liblfa.a -lsofa_c -lm

testobs: $(TESTOBS_OBJS) liblfa.a
	$(CC) -o testobs $(TESTOBS_OBJS) -L$(HOME)/lib64 liblfa.a -lsla -lm

testplan: $(TESTPLAN_OBJS) liblfa.a
	$(CC) -o testplan $(TESTPLAN_OBJS) liblfa.a -lncurses -lm

testrefro: $(TESTREFRO_OBJS) liblfa.a
	$(CC) -o testrefro $(TESTREFRO_OBJS) -L$(HOME)/lib64 liblfa.a -lsla -lm

testsun: $(TESTSUN_OBJS) liblfa.a
	$(CC) -o testsun $(TESTSUN_OBJS) -L$(HOME)/lib64 liblfa.a -lsla -lm

clean:
	rm -f $(OBJS) liblfa.a
	rm -f $(TESTJPL_OBJS) testjpl
	rm -f $(TESTOBS_OBJS) testobs
	rm -f $(TESTPLAN_OBJS) testplan
	rm -f $(TESTREFRO_OBJS) testrefro
	rm -f $(TESTSUN_OBJS) testsun
