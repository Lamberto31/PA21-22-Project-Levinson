#include <stdio.h>
#include <stdlib.h>

#define LOOP_COUNT = 1

double levinson(double*, double*, long);

int main(int argc, char const *argv[]) {
  double t[] = { 6, 4, 2, 1, 3, 5, 7 };
  double y[] = { 1, 2, 3,4};
  double n = 4;
  /*double t[] = { 2, 1, 3 };
  double y[] = { 1, 2};
  double n = 2;*/

  levinson(t, y, n);
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

  //Variabili necessarie per aggiornamento f e b contemporaneo
  double f_temp;
  double b_temp;

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
      b_temp = alpha_b * f[i] + beta_b * b[n-1-it+i];

      f[i] = f_temp;
      b[n-1-it+i] = b_temp;

      x[i] = x[i] + ((y[it] - e_x) * b[n-1-it+i]);
    }
  }
  for (int i = 0; i < n; i++) {
    fprintf(stdout, "x[%d]=%f\n", i, x[i]);
  }
}
