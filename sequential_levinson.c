#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_VALUE 10
#define LOOP_COUNT 10

double levinson(double*, double*, long);
void random_vector_generator(long, double*, int);

int main(int argc, char const *argv[]) {

  if(argc != 2){
 	  fprintf(stderr, "Usage: %s <n>\n", argv[0]);
    return -1;
  }
  long n = strtol(argv[1], NULL, 10);
  double *t;
  double *y;

  t = (double *) calloc(2*n-1, sizeof(double));
  y = (double *) calloc(n, sizeof(double));

  srand(time(NULL));
  while (!t[n-1]) {
    //TODO: controllare anche se tutti uguali???
    random_vector_generator(2*n-1, t, MAX_VALUE);
  }
  random_vector_generator(n, y, MAX_VALUE);

  levinson(t, y, n);
  free(t), t = NULL;
  free(y), y = NULL;
  return 0;
}

double levinson(double *t, double *y, long n) {

  //DICHIARAZIONE VARIABILI
  //Vettore avanti e indietro
  double *f;
  double *b;

  //Vettore delle incognite
  double *x;

  //Errori scalari per ogni estensione dei vettori f, b ed x
  double e_f;
  double e_b;
  double e_x;

  //Variabili necessarie per correzione errore ad ogni iterazione
  double d;
  double alpha_f;
  double beta_f;
  double alpha_b;
  double beta_b;

  //Variabile temporanea per poter aggiornare b usando f
  double f_temp;

  //Variabili per benchmark
  int iterations;
  clock_t elapsed;

  f = (double *) calloc(n, sizeof(double));
  b = (double *) calloc(n, sizeof(double));
  x = (double *) calloc(n, sizeof(double));

  elapsed = -clock();
  //TODO: inizializzare vettori dato che dopo la prima iterazione sono pieni
  for(iterations = 0; iterations < LOOP_COUNT; iterations++) {
    //CASO BASE
    f[0] = 1/t[n-1];
    b[n-1] = 1/t[n-1];
    x[0] = y[0]/t[n-1];

    for (int it = 1; it < n; it++) {

      e_f = 0;
      e_b = 0;
      e_x = 0;

      for (int i = 0; i < it; i++) {
        e_f = e_f + t[(it+1)-(i+1)+n-1] * f[i];
        e_b = e_b + t[-(i+1)+n-1] * b[n-it+i];
        e_x = e_x + t[(it+1)-(i+1)+n-1] * x[i];
      }

      d = 1 - (e_f * e_b);
      alpha_f = 1/d;
      beta_f = -e_f/d;
      alpha_b = -e_b/d;
      beta_b = 1/d;

      for (int i = 0; i < it+1; i++) {

        f_temp = alpha_f * f[i] + beta_f * b[n-1-it+i];
        b[n-1-it+i] = alpha_b * f[i] + beta_b * b[n-1-it+i];
        f[i] = f_temp;

        x[i] = x[i] + ((y[it] - e_x) * b[n-1-it+i]);
      }
    }
  }
  elapsed += clock();
  //TEST
  for (int i = 0; i < 2*n-1; i++) {
    fprintf(stdout, "t[%d] = %f\n", i, t[i]);
  }
  for (int i = 0; i < n; i++) {
    fprintf(stdout, "y[%d] = %f\n", i, y[i]);
  }
  for (int i = 0; i < n; i++) {
    fprintf(stdout, "x[%d] = %f\n", i, x[i]);
  }
  fprintf(stderr, "Tempo in media: %10.10lf Iterazioni: %d\n", ((double) elapsed / (double) iterations), iterations);
  //ENDTEST
  free(f), f = NULL;
  free(b), b = NULL;
  free(x), x = NULL;
}

void random_vector_generator(long n, double *v, int max) {
  for (long i = 0; i < n; i++) {
    v[i] = rand() % (max+1);
  }
}
