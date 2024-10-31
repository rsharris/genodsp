CFLAGS = -O3 -Wall -Wextra -Werror
LDLIBS = -lm

operators = sum clump percentile add multiply mask logical minmax morphology map opio variables

incFiles   = utilities.h genodsp_interface.h
opIncFiles = $(foreach op,${operators},${op}.h)

default: genodsp

genodsp: genodsp.o utilities.o $(foreach op,${operators},${op}.o)

%.o: %.c Makefile ${incFiles} ${opIncFiles}
	${CC} -c ${CFLAGS} $< -o $@

clean: cleano
	rm -f genodsp

cleano:
	rm -f *.o

install:
	mv genodsp $${HOME}/bin

archive:
	rm -f  genodsp.tar.gz
	mkdir  genodsp-distrib
	cp README              genodsp-distrib/
	cp genodsp.c           genodsp-distrib/
	cp genodsp_interface.h genodsp-distrib/
	cp Makefile            genodsp-distrib/
	cp add.c               genodsp-distrib/
	cp add.h               genodsp-distrib/
	cp clump.c             genodsp-distrib/
	cp clump.h             genodsp-distrib/
	cp logical.c           genodsp-distrib/
	cp logical.h           genodsp-distrib/
	cp map.c               genodsp-distrib/
	cp map.h               genodsp-distrib/
	cp mask.c              genodsp-distrib/
	cp mask.h              genodsp-distrib/
	cp minmax.c            genodsp-distrib/
	cp minmax.h            genodsp-distrib/
	cp morphology.c        genodsp-distrib/
	cp morphology.h        genodsp-distrib/
	cp multiply.c          genodsp-distrib/
	cp multiply.h          genodsp-distrib/
	cp opio.c              genodsp-distrib/
	cp opio.h              genodsp-distrib/
	cp percentile.c        genodsp-distrib/
	cp percentile.h        genodsp-distrib/
	cp sum.c               genodsp-distrib/
	cp sum.h               genodsp-distrib/
	cp utilities.c         genodsp-distrib/
	cp utilities.h         genodsp-distrib/
	cp variables.c         genodsp-distrib/
	cp variables.h         genodsp-distrib/
	rm -f genodsp-distrib/._*   # remove mac osx hidden files
	tar -cf - genodsp-distrib | gzip > genodsp.tar.gz
	rm -rf genodsp-distrib
