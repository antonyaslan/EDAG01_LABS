#include "simplex.h"

// Function to initialize the simplex struct
simplex_t* allocate_simplex(int m, int n) {
    simplex_t *s = (simplex_t *)malloc(sizeof(simplex_t));
    if (s == NULL) return NULL;

    s->m = m;
    s->n = n;

    // Allocate memory for variable arrays
    s->var = (int *)calloc((n + m + 1), sizeof(int));
    s->a = (double **)calloc(m, sizeof(double *));
    for (int i = 0; i < m; i++) {
        s->a[i] = (double *)calloc((n + 2), sizeof(double)); // Allocate extra column
    }
    s->b = (double *)calloc(m, sizeof(double));
    s->x = (double *)calloc((n + 1), sizeof(double));
    s->c = (double *)calloc(n, sizeof(double));

    s->y = 0.0; // Initialize y to 0

    return s;
}

// Function to free the simplex struct
void free_simplex(simplex_t *s) {
    if (s == NULL) return;

    free(s->var);
    for (int i = 0; i < s->m; i++) {
        free(s->a[i]);
    }
    free(s->a);
    free(s->b);
    free(s->x);
    free(s->c);
    free(s);
}

int init_simplex(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
    s->m = m;
    s->n = n;
    s->a = a;
    s->b = b;
    s->c = c;
    s->x = x;
    s->y = y;

    if(s->var == NULL) {
        printf("Failed to allocate memory for var");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < m + n + 1; i++) {
        s->var[i] = i;
    }

    int k = 0;
    for(int i = 1; i < m; i++) {
        if( b[i] < b[k]) {
            k = i;
        }
    }
    return k;
}

int select_nonbasic(simplex_t* s) {
    for (int i = 0; i < s->n; i++) {
        if (s->c[i] > E) {
            return i;
        }
    }
    return -1;
}

void pivot(simplex_t* s, int row, int col) {
    double** a = s->a;
    double* b = s->b;
    double* c = s->c;
    int m = s->m;
    int n = s->n;
    
    int t = s->var[col];
    s->var[col] = s->var[n + row];
    s->var[n + row] = t;
    s->y = s->y + c[col] * b[row] / a[row][col];

    for(int i = 0; i < n; i++) {
        if (i != col) {
            c[i] = c[i] - c[col] * a[row][i] / a[row][col];
        }
    }

    c[col] = -c[col] / a[row][col];

    for(int i = 0; i < m; i ++) {
        if (i != row) {
            for(int j = 0; j < n; j++) {
                if (j != col) {
                    a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
                }
            }
        }
    }

    for(int i = 0; i < m; i++) {
        if(i != row) {
            a[i][col] = -a[i][col] / a[row][col];
        }
    }

    for(int i = 0; i < n; i++) {
        if(i != col) {
            a[row][i] = a[row][i] / a[row][col];
        }
    }

    b[row] = b[row] / a[row][col];
    a[row][col] = 1.0 / a[row][col];
}

void prepare(simplex_t* s, int k) {
    int m = s->m;
    int n = s->n;

    for(int i = m + n; i > n; i--) {
        s->var[i] = s->var[i-1];
    }
    s->var[n] = m + n;
    n = n + 1;
    for(int i = 0; i < m; i++) {
        s->a[i][n-1] = -1;
    }
    s->x = (double*)calloc(m + n, sizeof(double));
    s->c = (double*)calloc(n, sizeof(double));

    s->c[n-1] = -1;
    s->n = n;
    pivot(s, k, n - 1);
}

int initial(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
    double w;
    int i,j,k;
    k = init(s, m, n, a, b, c, x, y, var);

    if(b[k] >= 0) {
        return 1; // feasible
    }
    prepare(s, k);
    n = s->n;
    s->y = xsimpelx(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);

    for(i=0; i < m + n; i++) {
        if(s->var[i] == m + n -1) {
            if (abs(s->x[i]) > E) {
                free(s->x);
                free(s->c);
                return 0; // infeasible
            }
        } else {
            break;
        }
    }

    if (i >= n) {
        j = 0;
        for(k = 0; k < n; k++) {
            if(abs(s->a[i-n][k]) > abs(s->a[i-n][j])) {
                j = k;
            }
        }
        pivot(s, i -n, j);
        i = j;
    }

    if (i < n - 1) {
        k = s->var[i];
        s->var[i] = s->var[n-1];
        s->var[n-1] = k;
        for(k =0; k < m; k = k +1) {
            w = s->a[k][n=1];
            s->a[k][n-1] = s->a[k][i];
            s->a[k][i] = w;
        }
    } else {
        
    }
    free(s->c);
    s->c = (double *)calloc(n, sizeof(double));
    s->c = c;
    s->y = y;
    
    for(k = n-1; k < n + m -1; k = k + 1) {
        s->var[k] = s->var[k+1];
    }
    n = s->n = s->n -1;

    double* t = (double *)calloc(n, sizeof(double));
    
    for(k = 0; k < n; k++) {
        for(j = 0; j < n; j++) {
            if(k = s->var[j]) {
                t[j] = t[j] + s->c[k];
                goto nextk;
            }
        }

        for(j = 0; j < m; j++) {
            if(s->var[n+j] = k) {
                break;
            }
        }

        s->y = s->y + s->c[k] * s->b[j];

        for(i = 0; i < n; i++) {
            t[i] = t[i] - s->c[k] * s->a[j][i];
        }
        nextk:;
    }

    for(i = 0; i < n; i++) {
        s->c[i] = t[i];
    }
    free(t);
    free(s->x);
    return 1;
}

// xsimplex Function
double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h) {
    simplex_t* s;
    s = allocate_simplex(m, n);
    int i, row, col;

    if (initial(&s, m, n, a, b, c, x, y, var) == 0) {
        free(s->var);
        return NAN;
    }
    col = select_nonbasic(&s);
    while(col >= 0) {
        row = -1;
        for(i = 0; i < m; i++) {
            if (a[i][col] > E && row < 0 || (b[i] / a[i][col] < b[row] / a[row][col])) {
                row = i;
            }
        }
        if (row < 0) {
            free(s->var);
            return INT_MAX; // unbounded
        }
        pivot(&s, row, col);
        col = select_nonbasic(&s);
    }
    if(h == 0) {
        for(i = 0; i < n; i++) {
            if (s->var[i] < n) {
                x[s->var[i]] = 0;
            }
        }
        
        for(i = 0; i < m; i++) {
            if (s->var[n + i] < n) {
                x[s->var[n+i]] = s->b[i];
            }
        }
        free(s->var);
    } else {
        for(i = 0; i < n; i++) {
            x[i] = 0;
        }
        for(i = n; i < n + m; i++) {
            x[i] = s->b[i-n];
        }
    }
    double ret_y = s->y;
    free_simplex(s);
    return ret_y;
}

// Run simplex algorithm
double simplex(int m, int n, double** a, double* b, double* c, double* x, double y) {
    return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}