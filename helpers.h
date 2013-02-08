#ifndef __HELPERS_H
#define __HELPERS_H

#include <stdlib.h>

// Generar un aleatorio en [0,1]
#define ALEATORIO ((double)rand()/RAND_MAX)


// Constantes usadas
#define LIBRE 0
#define OXIGENO -1
#define MONOXIDO 1
#define IMPUREZA 2


# define L 	120
# define N	L*L
# define Cal	50000
# define pasos	100000
# define yx	0.1


#define POS(i, j) (i + L*j)
#define DECR(indice) ((indice+L-1)%L)
#define INCR(indice) ((indice+1)%L)

// Matriz con las matrices de adyacencias.
typedef struct {
  int cell[N];
  int u[N];
  int d[N];
  int l[N];
  int r[N];
  int ul[N];
  int ur[N];
  int dl[N];
  int dr[N];
} MATRIZ;

#endif
