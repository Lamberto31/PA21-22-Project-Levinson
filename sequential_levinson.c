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


}
