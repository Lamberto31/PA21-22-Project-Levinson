#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

#define MAX_VALUE 10
#define LOOP_COUNT 1

void random_vector_generator(long, double*, int);
void vector_t_split(long, double*, double*, double*);
void parallel_levinson(double*, double*, double*, double*, double*);

int main(int argc, char *argv[]) {

  //VARIABLES
  //Process handling
  int id;
  int p;

  //Input
  long n;
  long t_size;
  double *t;
  double t_0;
  double *y;

  //Decomposition
  long v_size;
  long xres_size;

  //Custom datatype
  MPI_Datatype interleaved_vector;
  MPI_Aint start_address;
  MPI_Aint address;
  MPI_Aint lb;
  MPI_Aint extent;
  MPI_Datatype interleaved_vector_resized;

  //Vectors
  double *f;
  double *b;
  double *buf;
  double *x;
  double *x_res;

  //Benchmark
  double elapsed_time;
  double max_time;
  int iterations;
  int loop_count;

  //MPI Initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  //Input check
  if(argc != 2 && argc != 3) {
    if(!id) {
      fprintf (stderr, "Usage: %s <n> [<loop count>]\nWhere:\n<n> is the order of the system of linear equation\n<loop_count> is the number of iterations wanted for benchmarking\n", argv[0]);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  n = strtol(argv[1], NULL, 10);
  t_size = (2*n)-1;
  if (argc == 3)
    loop_count = strtol(argv[2], NULL, 10);
  else
    loop_count = LOOP_COUNT;

  //Input reading made by p0
  if(!id) {
    t = (double *) calloc(t_size, sizeof(double));
    if(!t) {
        fprintf(stderr, "Processor %d: Not enough memory\n", id);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    y = (double *) calloc(n, sizeof(double));
    if(!y) {
        fprintf(stderr, "Processor %d: Not enough memory\n", id);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    //TEST: Random generation
    srand(time(NULL));
    while(!t[n-1]) {
      //TODO: controllare anche se tutti uguali??? Capire cosa causa nan
      random_vector_generator(t_size, t, MAX_VALUE);
    }
    t_0 = t[n-1];
    random_vector_generator(n, y, MAX_VALUE);
    //ENDTEST
  }


  //Decomposition
  v_size = n/p + (n%p !=0);   //Ceil
  xres_size = v_size*p;

  //Custom datatype
  MPI_Type_vector(v_size, 1, p, MPI_DOUBLE, &interleaved_vector);
  MPI_Type_commit(&interleaved_vector);

  MPI_Get_address(&t[0], &start_address);
  lb = MPI_Aint_diff(start_address, start_address);
  MPI_Get_address(&t[1], &address);
  extent = MPI_Aint_diff(address, start_address);

  MPI_Type_create_resized(interleaved_vector, lb, extent, &interleaved_vector_resized);
  MPI_Type_commit(&interleaved_vector_resized);

  //Input distribution
  if(id) {
    t = (double *) calloc(t_size, sizeof(double));
    if(!t) {
      fprintf(stderr, "Processor %d: Not enough memory\n", id);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    y = (double *) calloc(n, sizeof(double));
    if(!y) {
      fprintf(stderr, "Processor %d: Not enough memory\n", id);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  MPI_Bcast(t, t_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(y, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //Vectors initialization
  f = (double *) calloc(v_size, sizeof(double));
  if(!f) {
    fprintf(stderr, "Processor %d: Not enough memory\n", id);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  b = (double *) calloc(v_size, sizeof(double));
  if(!b) {
    fprintf(stderr, "Processor %d: Not enough memory\n", id);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  buf = (double *) calloc(v_size, sizeof(double));
  if(!buf) {
    fprintf(stderr, "Processor %d: Not enough memory\n", id);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  x = (double *) calloc(v_size, sizeof(double));
  if(!x) {
    fprintf(stderr, "Processor %d: Not enough memory\n", id);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  if(!id) {
    x_res = (double *) calloc(xres_size, sizeof(double));
  if(!x_res) {
    fprintf(stderr, "Processor %d: Not enough memory\n", id);
    MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  //Start the timer
	MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time = -MPI_Wtime();

  //Begin the iterations
  for(iterations = 0; iterations < loop_count; iterations++) {

    //Vectors reset necessary due to repeated iterations
    memset(f, 0, v_size*sizeof(double));
    memset(b, 0, v_size*sizeof(double));
    memset(buf, 0, v_size*sizeof(double));
    memset(x, 0, v_size*sizeof(double));
    if(!id)
    memset(x_res, 0, xres_size*sizeof(double));

    //TODO ALGORITHM

  }

  //Stop the timer
	MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time += MPI_Wtime();

  //Memory release and finalize
  free(t), t = NULL;
  free(y), y = NULL;

  free(f), f = NULL;
  free(b), b = NULL;
  free(buf), buf = NULL;
  free(x), x = NULL;

  if(!id)
    free(x_res), x_res = NULL;

  MPI_Type_free(&interleaved_vector);
  MPI_Type_free(&interleaved_vector_resized);

  MPI_Finalize();

  return 0;
}

void random_vector_generator(long n, double *v, int max) {
  for (long i = 0; i < n; i++) {
    v[i] = rand() % (max+1);
    //v[i] = i+1;
  }

  //TEST
  if(n==7) {
    v[0] = 6;
    v[1] = 4;
    v[2] = 2;
    v[3] = 1;
    v[4] = 3;
    v[5] = 5;
    v[6] = 7;
  }
  if(n==4) {
    v[0] = 1;
    v[1] = 2;
    v[2] = 3;
    v[3] = 4;
  }
  //ENDTEST
}

void vector_t_split(long n, double *t, double *t_p, double *t_n) {
  memcpy(t_n, t, n*sizeof(double)); //In reverse order
  memcpy(t_p, &t[n], n*sizeof(double));
  //TODO reverse???
}

void parallel_levinson(double *t, double *y, double *f, double *b, double *x) {

}
