CXX = gcc
DEGUB_FLAGS = -ggdb
OPTIMIZE_FLAGS = -O3
OMP_FLAG = -fopenmp
POSTFLAGS = -lm

sommerfeld:
	${CXX} ${OPTIMIZE_FLAGS} ${OMP_FLAG} main_sommerfeld.c -o output/exec ${POSTFLAGS}

qpe:
	${CXX} ${OPTIMIZE_FLAGS} ${OMP_FLAG} main_qpe.c -o output/exec ${POSTFLAGS}

debug:
	${CXX} ${DEGUB_FLAGS} ${OMP_FLAG} main_sommerfeld.c -o output/exec ${POSTFLAGS}

clean:
	rm -r output/*
	mkdir output/err output/out

update: ../main_sommerfeld.c ../main_qpe.c
	cp ../main_sommerfeld.c .
	cp ../main_qpe.c .