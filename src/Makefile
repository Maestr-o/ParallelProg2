CC=mpicc
CFLAGS=-Wall -Wextra
BUILD_SEQ=../build/prog_seq
BUILD_PARALLEL=../build/prog_parallel
LDFLAGS=-fopenmp -lm

all: seq parallel

seq: main_seq.o
	@$(CC) $(CFLAGS) $(LDFLAGS) main_seq.o -o $(BUILD_SEQ)

parallel: main_parallel.o
	@$(CC) $(CFLAGS) $(LDFLAGS) main_parallel.o -o $(BUILD_PARALLEL)

main_seq.o: main.c
	@$(CC) $(CFLAGS) $(LDFLAGS) -c main.c -D seq -o main_seq.o

main_parallel.o: main.c
	@$(CC) $(CFLAGS) $(LDFLAGS) -c main.c -D parallel -o main_parallel.o

clean:
	@rm -f *.o

rebuild:
	clean all