snyder:
	mpicc -O3 -o par.out parallel_levinson.c -Wall

clean:
	rm par.out
