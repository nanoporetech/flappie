INCLUDE=-I..
LIBRARY=-lblas -lm
CFLAGS=-std=c11 -O3 -march=native

.PHONY: all
all: read_crp write_random_crp

read_crp: read_crp.o scrappie_util.o scrappie_matrix.o
	${CC} -o $@ $^ ${LIBRARY}

write_random_crp: write_random_crp.o scrappie_util.o scrappie_matrix.o
	${CC} -o $@ $^ ${LIBRARY}

.PHONY: scrappie_matrix.o
scrappie_matrix.o:../scrappie_matrix.o

%.o:%.c
	${CC} ${CFLAGS} ${INCLUDE} -c $<

.PHONY: clean
clean:
	rm -f read_crp.o scrappie_util.o scrappie_matrix.o write_random_crp.o
