CXX = gcc
FLAGS = -ggdb -O3 -fopenmp
POSTFLAGS = -lm

sommerfeld:
	${CXX} ${FLAGS} main_sommerfeld.c -o output/exec ${POSTFLAGS}

qpe:
	${CXX} ${FLAGS} main_qpe.c -o output/exec ${POSTFLAGS}

clean:
	rm -r output/*
	mkdir output/err output/out

update: ../main_sommerfeld.c ../main_qpe.c
	cp ../main_sommerfeld.c .
	cp ../main_qpe.c .