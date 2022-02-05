#include <stdio.h>
#include <stdlib.h>

#define LOOP_COUNT 1
#define N 4

double levinson(double*, double*, long);
void random_vector_generator(long, double*);

int main(int argc, char const *argv[]) {
  /*double t[] = { 6, 4, 2, 1, 3, 5, 7 };
  double y[] = { 1, 2, 3,4};
  double n = 4;*/
  /*double t[] = { 2, 1, 3 };
  double y[] = { 1, 2};
  double n = 2;*/

  double t[N];
  double y[N];

  srand(time(NULL));
  random_vector_generator(2*N-1, &t);
  random_vector_generator(N, &y);

  levinson(t, y, N);
  return 0;
}

double levinson(double *t, double *y, long n){

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

  f = (double *) calloc(n, sizeof(double));
  b = (double *) calloc(n, sizeof(double));
  x = (double *) calloc(n, sizeof(double));

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
  for (int i = 0; i < 2*N-1; i++) {
    fprintf(stdout, "t[%d] = %f\n", i, t[i]);
  }
  for (int i = 0; i < N; i++) {
    fprintf(stdout, "y[%d] = %f\n", i, y[i]);
  }
  for (int i = 0; i < n; i++) {
    fprintf(stdout, "x[%d]=%f\n", i, x[i]);
  }
}

void random_vector_generator(long n, double *v) {
  for (long i = 0; i < n; i++) {
    v[i] = rand() % 100;
  }
}
