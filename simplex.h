#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define E 1e-6
//Simplex Struct Definition

/**
 * @brief Represents a simplex stucture.
 * 
 * The struct stores the matrix (A), vectors (b, x, c), and the objective function value (y).
 * It also tracks the number of coefficients (m) and decision variables (n).
 */

typedef struct simplex_t {
    int m;          // Number of coefficients
    int n;          // Number of decision variables
    int *var;       // Array for variable indices 0..n-1 are nonbasic. size: n + m + 1
    double **a;     // Matrix A. size: m x (n + 2)
    double *b;      // Vector b. size: m
    double *x;      // Vector x. size: n + 1
    double *c;      // Vecotr c. size: n
    double y;       // Objective function value y
} simplex_t;

// Function Declarations

/**
 * @brief Allocates and initializes a simplex structure.
 * 
 * This function dynamically allocates memory for the simplex tableau, including
 * the matrix (A), vectors (b, x, c), and the variable index array (var). It also
 * ensures that the matrix (A) has an extra column for handling the additional variable.
 * 
 * @param int m Number of coefficients (rows in the matrix A).
 * @param int n Number of decision variables (columns in the matrix A without the extra column).
 * @return simplex_t Pointer to the initialized simplex structure, or NULL on failure.
 * 
 * @note The caller is responsible for freeing the allocated memory using free_simplex.
 */

simplex_t* allocate_simplex(int m, int n);

/**
 * @brief Frees the memory associated with a simplex structure.
 * 
 * This function deallocates all dynamically allocated memory in the simplex tableau,
 * including the matrix (A), vectors (b, x, c), and the variable index array (var).
 * 
 * @param simplex_t Pointer to the simplex structure to be freed.
 * 
 */

void free_simplex(simplex_t *s);

int init_simplex(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var);

int select_nonbasic(simplex_t* s);

void pivot(simplex_t* s, int row, int col);

void prepare(simplex_t* s, int k);

int initial(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var);

double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);

double simplex(int m, int n, double** a, double* b, double* c, double* x, double y);

#endif //SIMPLEX_H