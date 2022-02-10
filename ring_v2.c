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
  int *b;
  int b_size;
  int *buf;
  int ring_size;

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

  //Allocation of b e buf
  b_size = n/p;
  if(!id)
    fprintf(stdout, "b_size=%d\n", b_size);
  b = (int *) calloc(b_size, sizeof(int));
  if(!b){
        fprintf(stderr, "Processor %d: Not enough memory\n", id);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    buf = (int *) calloc(b_size, sizeof(int));
    if(!buf){
          fprintf(stderr, "Processor %d: Not enough memory\n", id);
          MPI_Abort(MPI_COMM_WORLD, -1);
      }

    for (int i = 0; i < b_size; i++) {
      b[i] = id + i*p;
      //b[i] = id;
    }
    /*
    fprintf(stdout, "IT = 0\n");
    for (int i = 0; i < b_size; i++) {
      fprintf(stdout, "id = %d\tb[%d] = %d\n", id, i, b[i]);
    }
    fprintf(stdout, "\n");
    */

    //Iterations
    for (int it = 1; it < n+1; it++) {

      if (id < it+1) {
        if (it+1 < p) {
          ring_size = it+1;
        } else {
          ring_size = p;
        }
        memcpy(buf, b, b_size*sizeof(int));
        if (id != 0)
          MPI_Recv(b, b_size, MPI_INT, id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Send(buf, b_size, MPI_INT, (id + 1) % ring_size, 0, MPI_COMM_WORLD);

        // Now process 0 can receive from the last process.
        if (id == 0)
          MPI_Recv(b, b_size, MPI_INT, ring_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      }
      if(id == p-1) {
        fprintf(stdout, "IT = %d\n", it);
        for (int i = 0; i < b_size; i++) {
          fprintf(stdout, "id = %d\tb[%d] = %d\n", id, i, b[i]);
        }
        fprintf(stdout, "\n");
      }

    }

  //Memory release and finalize
  free(b), b = NULL;
  free(buf), buf = NULL;
  MPI_Finalize();
  return 0;
}
