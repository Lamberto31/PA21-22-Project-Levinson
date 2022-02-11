#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
  int id;
  int p;

  int val;
  MPI_Status status;
  MPI_Request request;

  int n;
  int *x_res;
  int *x;
  int x_size;

  MPI_Datatype interleaved_vector;
  MPI_Datatype interleaved_vector_resized;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  //Input check
  if (argc != 2 && argc != 3) {
    if(!id) {
      fprintf (stderr, "Usage: %s <n>\nWhere:\n<n> is the order of the system of linear equation\n", argv[0]);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  n = strtol(argv[1], NULL, 10);

  //Allocation of x and x_res
  x_size = (n/p + (n%p !=0));
  if(!id) {
    fprintf(stdout, "x_size=%d\n", x_size);
    x_res = (int *) calloc(x_size*p, sizeof(int));
    if(!x_res){
          fprintf(stderr, "Processor %d: Not enough memory\n", id);
          MPI_Abort(MPI_COMM_WORLD, -1);
      }
  }
  x = (int *) calloc(x_size, sizeof(int));
  if(!x){
        fprintf(stderr, "Processor %d: Not enough memory\n", id);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    //Vector filling
    for (int i = 0; i < x_size; i++) {
      x[i] = id + i*p;
    }

    //Vector print
    for (int i = 0; i < x_size; i++) {
      fprintf(stdout, "id = %d\tx[%d] = %d\n", id, i, x[i]);
    }
    fprintf(stdout, "\n");

    //Custom datatype
    MPI_Type_vector(x_size, 1, p, MPI_INT, &interleaved_vector);
    MPI_Type_commit(&interleaved_vector);

    MPI_Type_create_resized(interleaved_vector, 0, sizeof(int), &interleaved_vector_resized);
    MPI_Type_commit(&interleaved_vector_resized);

    //Gather
    MPI_Gather(x, x_size, MPI_INT, x_res, 1, interleaved_vector_resized, 0, MPI_COMM_WORLD);

    //Result print
    if(!id) {
      for (int i = 0; i < n; i++) {
        fprintf(stdout, "x_res[%d] = %d\n", i, x_res[i]);
      }
      fprintf(stdout, "\n");
    }

  //Memory release and finalize
  free(x), x = NULL;
  if(!id)
    free(x_res), x_res = NULL;
  MPI_Type_free(&interleaved_vector);
  MPI_Type_free(&interleaved_vector_resized);
  MPI_Finalize();
  return 0;
}
