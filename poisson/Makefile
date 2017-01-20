CC = mpicc
BGPTASKNAME = bgp_run.jcf
NTHREADS = 4

main: main.c
	$(CC) -o main main.c -lm -std=c99 -fopenmp -O3

llsubmit: main
	llsubmit $(BGPTASKNAME)

runlocal: main
	mpirun -n $(NTHREADS) -host localhost main

clean:
	rm main
