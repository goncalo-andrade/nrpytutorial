CXX = gcc
DEGUB_FLAGS = -ggdb
OPTIMIZE_FLAGS = -O3
OMP_FLAG = -fopenmp
POSTFLAGS = -lm

sommerfeld:
	cp ../main_sommerfeld.c .
	${CXX} ${OPTIMIZE_FLAGS} ${OMP_FLAG} main_sommerfeld.c -o output/exec ${POSTFLAGS}

qpe:
	cp ../main_qpe.c .
	${CXX} ${OPTIMIZE_FLAGS} ${OMP_FLAG} main_qpe.c -o output/exec ${POSTFLAGS}

compile_only:
	${CXX} ${OPTIMIZE_FLAGS} ${OMP_FLAG} main_sommerfeld.c -o output/exec ${POSTFLAGS}

clean:
	rm -r output/*
	mkdir output/err output/out