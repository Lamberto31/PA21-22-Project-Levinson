#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

#define MAX_VALUE 10
#define LOOP_COUNT 1

void allocate_and_check(double**, long, int);
double random_input_generator(long, long, double*, double*);
void random_vector_generator(long, double*, int);
void vector_t_split(long, double*, double*, double*);
void create_resized_interleaved_vector_datatype(long, int, MPI_Datatype*);
void divide_work(long, int, int, long*, long*, int*, long*);
void exchange_vector(int, int, double*, long);
void parallel_levinson(int, int, long, double*, double*, long, double*, double*, double*, double []);
void print_toeplitz_matrix(long n, double*);
void print_result(long, long, double*, double*, double*, double, int, double[]);

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

  //Vectors
  double *f;
  double *b;
  double *x;
  double *x_res;

  //Benchmark
  double elapsed_time;
  double max_time;
  int iterations;
  int loop_count;
  double communication_time[2];
  double communication_max_time[2];

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
    allocate_and_check(&t, t_size, id);
    allocate_and_check(&y, n, id);
    t_0 = random_input_generator(n, t_size, t, y);
  }

  //Decomposition
  v_size = n/p + (n%p !=0);   //Ceil
  xres_size = v_size*p;

  //Custom datatype
  create_resized_interleaved_vector_datatype(v_size, p, &interleaved_vector);

  //Input distribution
  if(id) {
    allocate_and_check(&t, t_size, id);
    allocate_and_check(&y, n, id);
  }

  MPI_Bcast(t, t_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(y, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //Vectors initialization
  allocate_and_check(&f, v_size, id);
  allocate_and_check(&b, v_size, id);
  allocate_and_check(&x, v_size, id);
  if(!id)
    allocate_and_check(&x_res, xres_size, id);

  //Start the timer
	MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time = -MPI_Wtime();
  communication_time[0] = 0;
  communication_time[1] = 0;

  //Begin the iterations
  for(iterations = 0; iterations < loop_count; iterations++) {

    //Vectors reset necessary due to repeated iterations
    memset(f, 0, v_size*sizeof(double));
    memset(b, 0, v_size*sizeof(double));
    memset(x, 0, v_size*sizeof(double));
    if(!id)
    memset(x_res, 0, xres_size*sizeof(double));

    //EXECUTION OF ALGORITHM
    //Base case done just by process 0
    if(!id) {
      f[0] = 1/t_0;
      b[0] = 1/t_0;
      x[0] = y[0]/t_0;
    }

    //Call of parallel_levinson()
    parallel_levinson(id, p, n, t, y, v_size, f, b, x, communication_time);

  }

  //Stop the timer
	MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time += MPI_Wtime();

  //Gather
  MPI_Gather(x, v_size, MPI_DOUBLE, x_res, 1, interleaved_vector, 0, MPI_COMM_WORLD);

  //Reduction for max time
  MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&communication_time, &communication_max_time, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  //Result print
  if(!id)
    print_result(n, t_size, t, y, x_res, max_time, iterations, communication_time);

  //Memory release and finalize
  free(t), t = NULL;
  free(y), y = NULL;

  free(f), f = NULL;
  free(b), b = NULL;
  free(x), x = NULL;

  if(!id)
    free(x_res), x_res = NULL;

  MPI_Type_free(&interleaved_vector);

  MPI_Finalize();

  return 0;
}

void allocate_and_check(double **v, long v_size, int id) {
  *v = (double *) calloc(v_size, sizeof(double));
  if(!*v) {
    fprintf(stderr, "Processor %d: Not enough memory\n", id);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
}

double random_input_generator(long n, long t_size, double *t, double *y) {
  srand(time(NULL));
  while(!t[n-1]) {
    random_vector_generator(t_size, t, MAX_VALUE);
  }
  random_vector_generator(n, y, MAX_VALUE);
  return t[n-1];
}

void random_vector_generator(long n, double *v, int max) {
  for (long i = 0; i < n; i++) {
    v[i] = rand() % max + 1;
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
    //TEST
    if(n==19) {
      v[0] = 7;
      v[1] = 2;
      v[2] = 5;
      v[3] = 9;
      v[4] = 7;
      v[5] = 4;
      v[6] = 4;
      v[7] = 4;
      v[8] = 7;
      v[9] = 4;
      v[10] = 7;
      v[11] = 5;
      v[12] = 4;
      v[13] = 8;
      v[14] = 8;
      v[15] = 9;
      v[16] = 9;
      v[17] = 7;
      v[18] = 5;
    }
    if(n==10) {
      v[0] = 4;
      v[1] = 2;
      v[2] = 10;
      v[3] = 5;
      v[4] = 9;
      v[5] = 7;
      v[6] = 10;
      v[7] = 4;
      v[8] = 1;
      v[9] = 9;
    }
    //ENDTEST
}

void create_resized_interleaved_vector_datatype(long n, int stride, MPI_Datatype *resized_interleaved_vector_datatype) {
  MPI_Datatype type_vector;

  double vector[2];
  MPI_Aint start_address;
  MPI_Aint address;
  MPI_Aint lb;
  MPI_Aint extent;

  MPI_Type_vector(n, 1, stride, MPI_DOUBLE, &type_vector);
  MPI_Type_commit(&type_vector);

  MPI_Get_address(&vector[0], &start_address);
  lb = MPI_Aint_diff(start_address, start_address);
  MPI_Get_address(&vector[1], &address);
  extent = MPI_Aint_diff(address, start_address);

  MPI_Type_create_resized(type_vector, lb, extent, resized_interleaved_vector_datatype);
  MPI_Type_commit(resized_interleaved_vector_datatype);

  MPI_Type_free(&type_vector);
}


void divide_work(long it, int id, int p, long *ops_errors, long *ops_update, int *ring_size, long *els_to_exchange) {

  *ops_errors = it/p + (it % p > id);
  *ops_update = (it+1)/p + ((it+1) % p > id);
  if (it+1 < p) {
    *ring_size = it+1;
  } else {
    *ring_size = p;
  }
  *els_to_exchange = it/p + (it%p !=0);
}

void exchange_vector(int ring_size, int id, double *v, long v_size) {
  if (ring_size > 1) {
    double *buf;
    MPI_Request req;

    allocate_and_check(&buf, v_size, id);

    memcpy(buf, v, v_size*sizeof(double));
    if (!id) {
      MPI_Irecv(v, v_size, MPI_DOUBLE, ring_size - 1, 0, MPI_COMM_WORLD, &req);
    } else {
      MPI_Irecv(v, v_size, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD, &req);
    }
    MPI_Send(buf, v_size, MPI_DOUBLE, (id + 1) % ring_size, 0, MPI_COMM_WORLD);
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    free(buf), buf = NULL;
  }
}

void parallel_levinson(int id, int p, long n, double *t, double *y, long v_size, double *f, double *b, double *x, double comm_time[]) {

  //Work division
  long ops_errors;
  long ops_update;
  int ring_size;
  long els_to_exchange;

  //Errors (0 is e_f, 1 is e_b, 2 is e_x)
  double errors[3];
  double global_errors[3];

  //Correctors
  double d;
  double alpha_f;
  double beta_f;
  double alpha_b;
  double beta_b;
  double beta_x;

  //Vectors Update
  double f_temp;

  for (long it = 1; it < n; it++) {
    divide_work(it, id, p, &ops_errors, &ops_update, &ring_size, &els_to_exchange);

    //Errors initialization and computation
    memset(errors, 0, 3*sizeof(double));

    for (long i = 0; i < ops_errors; i++) {
      errors[0] += t[it-id-i*p+n-1] * f[i];
      errors[1] += t[-id-1-i*p+n-1] * b[ops_errors-1-i];
      errors[2] += t[it-id-i*p+n-1] * x[i];
    }

    //Reduction
  	MPI_Barrier(MPI_COMM_WORLD);
    comm_time[0] += -MPI_Wtime();
    MPI_Allreduce(&errors, &global_errors, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    comm_time[0] += MPI_Wtime();

    //Correctors computation
    d = 1 - (global_errors[0] * global_errors[1]);
    alpha_f = 1/d;
    beta_f = -global_errors[0]/d;
    alpha_b = -global_errors[1]/d;
    beta_b = 1/d;
    beta_x = y[it] - global_errors[2];


    //Vector b exchange
    MPI_Barrier(MPI_COMM_WORLD);
    comm_time[1] += -MPI_Wtime();
    if (ops_update)
      exchange_vector(ring_size, id, b, els_to_exchange);
    MPI_Barrier(MPI_COMM_WORLD);
    comm_time[1] += MPI_Wtime();

    //Vectors f,b,x update{
      for (long i = 0; i < ops_update; i++) {
        f_temp = alpha_f * f[i] + beta_f * b[ops_update-1-i];
        b[ops_update-1-i] = alpha_b * f[i] + beta_b * b[ops_update-1-i];
        f[i] = f_temp;
        x[i] = x[i] + (beta_x * b[ops_update-1-i]);
      }
  }
}

void print_toeplitz_matrix(long n, double *t) {
  long half = (n-1)/2;
  if (n > 5) {
    for (long i = 0; i < (2*n)-1; i++) {
      fprintf(stdout, "t[%ld] = %10.10lf\n", i, t[i]);
    }
  } else {
    for (long i = 0; i < n; i++) {
      if (i == half)
        fprintf(stdout, "T =");
      fprintf(stdout, "\t[\t");
      for (long j = 0; j < n; j++) {
        fprintf(stdout, "%f\t", t[i-j+n-1]);
      }
      fprintf(stdout, "]\n");
    }
  }
}

void print_result(long n, long t_size, double *t, double *y, double *x_res, double time, int iterations, double comm_time[]) {
  double average_time;
  double average_reduction_time;
  double average_exchange_vector_time;
  double average_communication_time;
  double average_fraction_communication_total_time;

  if (n<50) {
    print_toeplitz_matrix(n, t);
    fprintf(stdout, "\n");
    for(long i = 0; i < n; i++) {
      fprintf(stdout, "y[%ld] = %10.10lf\n", i, y[i]);
    }
    fprintf(stdout, "\n");
    for(long i = 0; i < n; i++) {
      fprintf(stdout, "x_res[%ld] = %10.10lf\n", i, x_res[i]);
    }
    fprintf(stdout, "\n");

    fprintf(stderr, "t = [");
    for (long i = 0; i < (2*n)-1; i++) {
      fprintf(stderr, "%0.0lf", t[i]);
      if (i != (2*n)-2) {
        fprintf(stderr, ",");
      }
    }
    fprintf(stderr, "]\n");

    fprintf(stderr, "y = [");
    for (long i = 0; i < n; i++) {
      fprintf(stderr, "%0.0lf", y[i]);
      if (i != n-1) {
        fprintf(stderr, ",");
      }
    }
    fprintf(stderr, "]\n");
  }

  average_time = (double) time / (double) iterations;
  average_reduction_time = (double) comm_time[0] / (double) iterations;
  average_exchange_vector_time = (double) comm_time[1] / (double) iterations;
  average_communication_time = average_reduction_time + average_exchange_vector_time;
  average_fraction_communication_total_time = average_communication_time / average_time *100;

  fprintf(stderr, "Average total time: \t\t\t%10.10lf\tIterazioni: %d\n", average_time, iterations);
  fprintf(stderr, "Average all_reduction time: \t\t%10.10lf\tIterazioni: %d\n", average_reduction_time, iterations);
  fprintf(stderr, "Average exchange vector time: \t\t%10.10lf\tIterazioni: %d\n", average_exchange_vector_time, iterations);
  fprintf(stderr, "Average communication time: \t\t%10.10lf\tIterazioni: %d\n", average_communication_time, iterations);
  fprintf(stderr, "Fraction of communication time: \t%0.2f%%\t\tIterazioni: %d\n", average_fraction_communication_total_time, iterations);

}
