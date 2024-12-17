#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double** make_matrix(int m, int n) {
    double**    matrix;
    int         i;

    matrix = calloc(m,sizeof(double*));
    for (i = 0; i < m; i += 1) {
        matrix[i] = calloc(n , sizeof(double));
    }
    return matrix;
}

double* make_coeff_vector(int n) {
    double* coeff_vector = calloc(n , sizeof(double));
    return coeff_vector;
}

double* make_constants(int m) {
    double* constants = calloc(m, sizeof(double));
    return constants;
}

void print_coeff_equation(double* coeff_vector, int n) {
    printf("Coefficient Equation:\n");
    printf(" max z = ");
    for (int i = 0; i < n; i++) {
        printf("%10.3lf x%d ", coeff_vector[i], i);
        if (i < n -1 ) {
            printf("+");
        }
    }
    printf("\n");
}

void print_dec_variable_vector(double* dec_variable_vector, int m) {
    printf("Decision Variable Vector:\n");
    for (int i = 0; i < m; i++) {
        printf("%10.3lf ", dec_variable_vector[i]);
    }
    printf("\n");
}

void print_inequalities(double** matrix, double* dec_variable_vector, int m, int n) {
    printf("Inequalities:\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.3lf x%d ", matrix[i][j], j);
            if (j < n - 1) {
                printf("+");
            } else {
                printf("\u2264");
            }
        }
        printf(" %10.3lf", dec_variable_vector[i]);
        printf("\n");
    }
}

int main() {
    int m, n;
    double coeff, matrix_value, constant;
    //printf("Enter the size of matrix and  vectors: \n");
    scanf("%d %d", &m, &n);

    double** matrix = make_matrix(m,n);

    double* coefficients = make_coeff_vector(n);

    double* constants = make_constants(m);

    for(int j = 0; j < n; j+=1) {
        scanf("%lf", &coeff);
        coefficients[j] = coeff;
    }

    for(int i = 0; i < m; i+=1) {
        for(int j = 0; j< n; j+=1) {
            scanf("%lf", &matrix_value);
            matrix[i][j] = matrix_value;
        }
    }

    for(int i = 0; i < m; i+=1) {
        scanf("%lf", &constant);
        constants[i] = constant;
    }

    // Print the matrix and vectors
    print_coeff_equation(coefficients, n);
    print_inequalities(matrix, constants, m, n);

    // Free allocated memory
    for (int i = 0; i < m; i++) {
        free(matrix[i]);
    }
    free(matrix);
    free(coefficients);
    free(constants);

    return 0;
}
