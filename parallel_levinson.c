#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

#define MAX_VALUE 100
#define LOOP_COUNT 1000

double levinson(double*, double*, long);
void random_vector_generator(long, double*, int);

int main(int argc, char *argv[]) {

  //VARIABLES DECLARATIONS
  //Process handling
  int id;
  int p;

  //Input
  long n;
  double *t;
  //TODO Capire come gestire y uguale per tutti
  //double *y;
  //TEST
  /*
  long n = 4;
  double t[] = { 6, 4, 2, 1, 3, 5, 7 };
  double y[] = { 1, 2, 3,4};
  */
  //ENDTEST

  //Decomposition
  long vectors_size;
  long t_size;

  //Vectors
  double *f;
  double *b;
  double *x;

  //Iteration idx
  long it;

  //Errors for each extension of vectors f, b, x
  double e_f;
  double e_b;
  double e_x;

  //Correctors
  double d;
  double alpha_f;
  double beta_f;
  double alpha_b;
  double beta_b;

  //Temp variable for safe update of b and f
  double f_temp;

  //Benchmark
  double elapsed_time;
  double max_time;
  int iterations;
  int loop_count;

  //PRELIMINARY OPERATIONS
  //MPI Initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  //Input check
  if (argc != 2 && argc != 3) {
    if(!id) {
      fprintf (stderr, "Usage: %s <n> [<loop count>]\nWhere:\n<n> is the order of the system of linear equation\n<loop_count> is the number of iterations wanted for benchmarking\n", argv[0]);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  n = strtol(argv[1], NULL, 10);
  if (argc == 3)
    loop_count = strtol(argv[2], NULL, 10);
  else
    loop_count = LOOP_COUNT;

  //Decomposition
  vectors_size=n/p;
  t_size=2*vectors_size;

  //Input
  //TEST
  t = (double *) calloc(t_size, sizeof(double));
  //TODO Capire come gestire y uguale per tutti
  //y = (double *) calloc(n, sizeof(double));
  double y[] = { 1, 2, 3,4};

  srand(time(NULL));
  //TODO: capire come controllare t0
  /*while (!t[n-1]) {
    //TODO: controllare anche se tutti uguali??? Capire cosa causa nan
    random_vector_generator(t_size, t, MAX_VALUE);
  }*/
  random_vector_generator(t_size, t, MAX_VALUE);
  //random_vector_generator(n, y, MAX_VALUE);
  //ENDTEST

  //Vectors initialization
  f = (double *) calloc(vectors_size, sizeof(double));
  if(!f){
    fprintf(stderr, "Processor %d: Not enough memory\n", id);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  b = (double *) calloc(vectors_size, sizeof(double));
  if(!b){
    fprintf(stderr, "Processor %d: Not enough memory\n", id);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  if(!id)
    x = (double *) calloc(n, sizeof(double));
  else
    x = (double *) calloc(vectors_size, sizeof(double));
  if(!x){
    fprintf(stderr, "Processor %d: Not enough memory\n", id);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  //Start the timer
	MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time = -MPI_Wtime();

  //Begin the iterations
  for(iterations = 0; iterations < LOOP_COUNT; iterations++) {

    //Vectors reset necessary due to repeated iterations
    memset(f, 0, n*sizeof(double));
    memset(b, 0, n*sizeof(double));
    memset(x, 0, n*sizeof(double));

    //Base case done just by process 0
    if(!id) {
      f[0] = 1/t[n-1];
      b[0] = 1/t[n-1];
      x[0] = y[0]/t[n-1];
    }

    //Begin of the algorithm iterations
    for (it = 1; it < n; it++) {

      //Errors initialization and computation
      e_f = 0;
      e_b = 0;
      e_x = 0;
      for (long i = p; i < it; i+p) {
        e_f = e_f + t[(it+1)-(i+1)+n-1] * f[i];
        e_b = e_b + t[(i+1)-(it+1)+n-1] * b[i];
        e_x = e_x + t[(it+1)-(i+1)+n-1] * x[i];
      }

      //TODO ALL-REDUCE
      MPI_Barrier(MPI_COMM_WORLD);

      //Correctors computation (check to avoid useless work)
      if (id > it+1) {
        d = 1 - (e_f * e_b);
        alpha_f = 1/d;
        beta_f = -e_f/d;
        alpha_b = -e_b/d;
        beta_b = 1/d;
      }

      //Vectors update
      for (long i = p; i < it+1; i+p) {
        f_temp = alpha_f * f[i] + beta_f * b[it-i];
        b[it-i] = alpha_b * f[i] + beta_b * b[it-i];
        f[i] = f_temp;
        x[i] = x[i] + ((y[it] - e_x) * b[it-i]);
      }
    }
    //TODO GATHER
    MPI_Barrier(MPI_COMM_WORLD);
  }

  //Stop the timer
	MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time = -MPI_Wtime();

  //Result print
  //TEST
  if(!id){
    for (int i = 0; i < 2*n-1; i++) {
      fprintf(stdout, "t[%d] = %10.10lf\n", i, t[i]);
    }
    for (int i = 0; i < n; i++) {
      fprintf(stdout, "y[%d] = %10.10lf\n", i, y[i]);
    }
    for (int i = 0; i < n; i++) {
      fprintf(stdout, "x[%d] = %10.10lf\n", i, x[i]);
    }
    //TODO calcolo max_time con reduce;
    fprintf(stderr, "Tempo medio: %10.10lf Iterazioni: %d\n", ((double) elapsed_time / (double) iterations), iterations);
  }
  //ENDTEST

  //Memory release and finalize
  free(f), f = NULL;
  free(b), b = NULL;
  free(x), x = NULL;

  free(t), t = NULL;
  //free(y), y = NULL;

  MPI_Finalize();
  return 0;
}

void random_vector_generator(long n, double *v, int max) {
  for (long i = 0; i < n; i++) {
    //v[i] = rand() % (max+1);
    v[i] = i+1;
  }
}
