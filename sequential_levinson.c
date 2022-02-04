#include <stdio.h>
#include <stdlib.h>

#define LOOP_COUNT = 1

double levinson(double*, double*, long);

int main(int argc, char const *argv[]) {
  return 0;
}

double levinson(double *t, double *y, long n){

  //DICHIARAZIONE VARIABILI
  //Vettore avanti e indietro
  double f;
  double b;

  //Vettore delle incognite
  double x;

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


  f = (double *) calloc(n, sizeof(double));
  b = (double *) calloc(n, sizeof(double));
  x = (double *) calloc(n, sizeof(double));

  //CASO BASE
  f[0] = 1/t[n-1];
  b[0] = 1/t[n-1];
  x[0] = y[0]/t[n-1];

  for (int it = 1; i < n; i++) {
    f[it] = 0;
    b[it] = 0;
    x[it] = 0;

    e_f = 0;
    e_b = 0;
    e_x = 0;

    for (int i = 0; i < it; i++) {
      e_f = e_f + t[it-i+it] * f[i];
      e_b = e_b + t[-i+it] * b[i];
      e_x = e_x + t[it-i+it] * f[i];
    }
  }

}
