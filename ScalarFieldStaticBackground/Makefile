CXX = gcc
DEGUB_FLAGS = -ggdb
OPTIMIZE_FLAGS = -O3
OMP_FLAG = -fopenmp
POSTFLAGS = -lm

sommerfeld:
	${CXX} ${OPTIMIZE_FLAGS} ${OMP_FLAG} Sommerfeld_Playground/main_sommerfeld.c -o Sommerfeld_Playground/output/exec ${POSTFLAGS}

qpe:
	${CXX} ${OPTIMIZE_FLAGS} ${OMP_FLAG} QPE_Playground/main_qpe.c -o QPE_Playground/output/exec ${POSTFLAGS}

# debug:
# 	${CXX} ${DEGUB_FLAGS} ${OMP_FLAG} main_sommerfeld.c -o output/exec ${POSTFLAGS}

clean:
	rm -r */output/*
	mkdir Sommerfeld_Playground/output/err Sommerfeld_Playground/output/out
	mkdir QPE_Playground/output/err QPE_Playground/output/out

update: ../main_sommerfeld.c ../main_qpe.c
	cp main_sommerfeld.c Sommerfeld_Playground/
	cp main_qpe.c QPE_Playground/