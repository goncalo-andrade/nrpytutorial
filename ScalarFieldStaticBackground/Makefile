CXX = gcc
DEGUB_FLAGS = -ggdb
OPTIMIZE_FLAGS = -O3
OMP_FLAG = -fopenmp
POSTFLAGS = -lm

sommerfeld:
	${CXX} ${OPTIMIZE_FLAGS} ${OMP_FLAG} SommerfeldPlayground/main_sommerfeld.c -o SommerfeldPlayground/output/exec ${POSTFLAGS}

qpe:
	${CXX} ${OPTIMIZE_FLAGS} ${OMP_FLAG} QPEPlayground/main_qpe.c -o QPEPlayground/output/exec ${POSTFLAGS}

# debug:
# 	${CXX} ${DEGUB_FLAGS} ${OMP_FLAG} main_sommerfeld.c -o output/exec ${POSTFLAGS}

clean:
	rm -r */output/*
	mkdir SommerfeldPlayground/output/err SommerfeldPlayground/output/out
	mkdir QPEPlayground/output/err QPEPlayground/output/out

update: main_sommerfeld.c main_qpe.c
	cp main_sommerfeld.c SommerfeldPlayground/
	cp main_qpe.c QPEPlayground/