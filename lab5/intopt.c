#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#define E 1e-6

typedef struct simplex_t {
    int m;          // Number of constraints
    int n;          // Number of decision variables
    int *var;       // Array for variable indices 0..n-1 are nonbasic. size: n + m + 1
    double **a;     // Matrix A. size: m x (n + 2)
    double *b;      // Vector b. size: m
    double *x;      // Vector x. size: n + 1
    double *c;      // Vecotr c. size: n
    double y;       // Objective function value y
} simplex_t;

typedef struct node_t {
    int m;          // Number of constraints
    int n;          // Number of decision variables
    int k;          // Parent branches on xk
    int h;          // Branch on xh
    double xh;      // xh
    double ak;      // Parent ak
    double bk;      // Parent bk
    double *min;    // Lower bounds. size: n
    double *max;    // Upper bounds. size: n
    double **a;     // Matrix A. size: m x (n + 2)
    double *b;      // Vector b. size: m
    double *x;      // Vector x. size: n + 1
    double *c;      // Vecotr c. size: n
    double z;       // z
} node_t;

typedef struct set_t {
    node_t *data;
    struct set_t *next;

} set_t;

double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);

void print_simplex(simplex_t *s);

int init(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
    int i,k;

    s->m = m;
    s->n = n;
    s->a = a;
    s->b = b;
    s->c = c;
    s->x = x;
    s->y = y;
    s->var = var;

    if(s->var == NULL) {
        s->var = (int*) calloc(n + m + 1, sizeof(int)); 
        for (i = 0; i < m + n; i++) {
            s->var[i] = i;
        }
    }

    for(k=0, i = 1; i < m; i++) {
        if( b[i] < b[k]) {
            k = i;
        }
    }
    return k;
}

int select_nonbasic(simplex_t* s) {
    int i;
    for (i = 0; i < s->n; i++) {
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
    int i,j,t;
    t = s->var[col];
    s->var[col] = s->var[n + row];
    s->var[n + row] = t;
    s->y = s->y + c[col] * b[row] / a[row][col];

    for(i = 0; i < n; i++) {
        if (i != col) {
            c[i] = c[i] - c[col] * a[row][i] / a[row][col];
        }
    }

    c[col] = -c[col] / a[row][col];

    for(i = 0; i < m; i++) {
        if(i != row) {
            b[i] = b[i] - a[i][col] * b[row] / a[row][col];
        }
    }
    for(i = 0; i < m; i++) {
        if(i != row) {
            for(j = 0; j < n; j++) {
                if(j != col) {
                    a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
                }
            }
        }
    }

    for(i = 0; i < m; i++) {
        if(i != row) {
            a[i][col] = -a[i][col] / a[row][col];
        }
    }

    for(i = 0; i < n; i++) {
        if(i != col) {
            a[row][i] = a[row][i] / a[row][col];
        }
    }

    b[row] = b[row] / a[row][col];
    a[row][col] = 1.0 / a[row][col];
    //print_simplex(s);
}

void prepare(simplex_t* s, int k) {
    int m = s->m;
    int n = s->n;
    int i;
    for(i = m + n; i > n; i--) {
        s->var[i] = s->var[i-1];
    }
    s->var[n] = m + n;
    n = n + 1;
    for(i = 0; i < m; i++) {
        s->a[i][n-1] = -1;
    }
    s->x = (double*)calloc(m + n, sizeof(double));
    s->c = (double*)calloc(n, sizeof(double));
    s->c[n-1] = -1;
    s->n = n;
    pivot(s, k, n - 1);
}

int initial(simplex_t* s, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
    int i, j, k=0;
    double w;

    k = init(s, m, n, a, b, c, x, y, var);

    if (b[k] >= 0) {
      return 1; //feasible.
    }

    prepare(s, k);
    n = s->n;
    s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);

    for (i = 0; i < m + n; i++) {
        if (s->var[i] == m + n - 1) {
            if (fabs(s->x[i]) > E) {
                free(s->x);
                free(s->c);
                s->x = NULL;
                s->c = NULL;
                return 0; //infeasible.
            } else {
                break;
            }
        }
    }

    if (i >= n) {
        // x_{n+m} is basic. find good nonbasic.
        for (j=0, k=0; k < n; k++) {
            if (fabs(s->a[i - n][k]) > fabs(s->a[i - n][j])) {
                j = k;
            }
        }
        pivot(s, i - n, j);
        i = j;
    }

    if (i < n - 1) {
        //x_{n+m} is nonbasic and not last. swap columns i and n-1.
        k = s->var[i];
        s->var[i] = s->var[n - 1];
        s->var[n - 1] = k;
        for (k = 0; k < m; k++) {
            w = s->a[k][n - 1];
            s->a[k][n - 1] = s->a[k][i];
            s->a[k][i] = w;
        }
    } else {
        //x_{n+m} is nonbasic and last. forget it.
    }

    free(s->c);
    s->c = c;
    s->y = y;

    for (k = n - 1; k < n + m - 1; k++) {
        s->var[k] = s->var[k + 1];
    }

    n = s->n = s->n - 1;
    double *t = calloc(n, sizeof(double *));

    for (k = 0; k < n; k++) {
        for (j = 0; j < n; j++) {
           if (k == s->var[j]) {
               //x_k is nonbasic. add c_k.
               t[j] = t[j] + s->c[k];
               goto next_k;
           }
        }

        for (j = 0; j < m; j++) {
           if (s->var[n + j] == k) {
               //x_k is at row j.
               break;
           }
        }

        s->y = s->y + s->c[k] * s->b[j];

        for (i = 0; i < n; i++) {
            t[i] = t[i] - s->c[k] * s->a[j][i];
        }
        next_k:;
    }

    for (i = 0; i < n; i++) {
        s->c[i] = t[i];
    }
    free(t);
    free(s->x);
    t = NULL;
    s->x = NULL;
    return 1;
}

// xsimplex Function
double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h) {
    simplex_t s;
    int i, row, col;

    if (!initial(&s, m, n, a, b, c, x, y, var)) {
        free(s.var);
        s.var = NULL;
        return NAN;
    }
    while((col=select_nonbasic(&s)) >= 0) {
        row = -1;
        for(i = 0; i < m; i++) {
            if (a[i][col] > E && (row < 0 || b[i] / a[i][col] < b[row] / a[row][col])) {
                row = i;
            }
        }
        if (row < 0) {
            free(s.var);
            s.var = NULL;
            return INFINITY; // unbounded
        }
        pivot(&s, row, col);
    }

    if(h == 0) {
        for(i = 0; i < n; i++) {
            if (s.var[i] < n) {
                x[s.var[i]] = 0;
            }
        }
        
        for(i = 0; i < m; i++) {
            if (s.var[n + i] < n) {
                x[s.var[n+i]] = s.b[i];
            }
        }
        free(s.var);
        s.var = NULL;
    } else {
        for(i = 0; i < n; i++) {
            x[i] = 0;
        }
        for(i = n; i < n + m; i++) {
            x[i] = s.b[i-n];
        }
    }
    return s.y;
}

// Run simplex algorithm
double simplex(int m, int n, double** a, double* b, double* c, double* x, double y) {
    return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}

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
    double* coeff_vector = calloc(n + 1, sizeof(double));
    return coeff_vector;
}

double* make_constants(int m) {
    double* constants = calloc(m , sizeof(double));
    return constants;
}

double* make_dec_variables(int n) {
    double* dec_variable_vector = calloc(n, sizeof(double));
    return dec_variable_vector;
}


node_t * initial_node(int m, int n, double **a, double *b, double *c) {
    int i;
    node_t *p = (node_t *) calloc(1, sizeof(node_t));
    p->a = make_matrix(m+1, n+1);
    p->b = make_constants(m+1);
    p->c = make_coeff_vector(n);
    p->x = make_dec_variables(n+1);
    p->min = make_dec_variables(n);
    p->max = make_dec_variables(n);

    p->m = m;
    p->n = n;
    for (i = 0; i < m; i++) {
        memcpy(p->a[i], a[i], n * sizeof(double));
    }
    memcpy(p->b, b, m * sizeof(double));
    memcpy(p->c, c, n * sizeof(double));

    for (i = 0; i < n; i++) {
        p->min[i] = -INFINITY;
        p->max[i] = INFINITY;
    }

    return p;
}

node_t * extend(node_t *p, int m, int n, double **a, double *b, double *c, int k, double ak, double bk) {
    node_t *q = malloc(sizeof(*q));
    int i,j;
    q->k = k;
    q->ak = ak;
    q->bk = bk;

    if(ak > 0 && isfinite(p->max[k])) {
        q->m = p->m;
    } else if(ak < 0 && (p->min[k] > 0)) {
        q->m = p->m;
    } else {
        q->m = p->m + 1;
    }

    q->n = p->n;

    q->h = -1;

    q->a = make_matrix(q->m+1, q->n+1); // note normally q->m > m
    q->b = make_constants(q->m+1);
    q->c = make_coeff_vector(q->n);
    q->x = make_dec_variables(q->n+1);
    q->min = make_dec_variables(n);
    q->max = make_dec_variables(n);

    memcpy(q->min, p->min, n * sizeof(double));
    memcpy(q->max, p->max, n * sizeof(double));
    for (i = 0; i < m; i++) {
        memcpy(q->a[i], a[i], n * sizeof(double));
    }
    memcpy(q->b, b, m * sizeof(double));
    memcpy(q->c, c, n * sizeof(double));


    if (ak > 0) {
        if((q->max[k] == INFINITY) || bk < q->max[k]) {
            q->max[k] = bk;
        }
    } else if ((q->min[k] == -INFINITY) || (-bk > q->min[k])) {
        q->min[k] = -bk;
    }

    for(i=m, j=0; j<n; j++) {
        if(q->min[j] > -INFINITY) {
            q->a[i][j] = -1;
            q->b[i] = -q->min[j];
            i++;
        }
        if(q->max[j] < INFINITY) {
            q->a[i][j] = 1;
            q->b[i] = q->max[j];
            i++;
        }
    }
    return q;
}

int is_integer(double *xp) {
    double x = *xp;
    double r = lround(x); // ISO C lround
    if(fabs(r-x) < E) {
        *xp = r;
        return 1;
    } else {
        return 0;
    }
}

int integer(node_t *p) {
    int i;
    for(i=0; i<p->n; i++) {
        if(!is_integer(&(p->x[i]))) {
            return 0;
        }
    }
    return 1;
}

void free_node(node_t **p) {
    if (*p == NULL)
        return;
    if ((*p)->a != NULL) {
        for (int i = 0; i < (*p)->m + 1; i++) {
        free((*p)->a[i]);
        (*p)->a[i] = NULL;
        }
        free((*p)->a);
        (*p)->a = NULL;
    }
    if ((*p)->b != NULL) {
        free((*p)->b);
        (*p)->b = NULL;
    }
    if ((*p)->c != NULL) {
        free((*p)->c);
        (*p)->c = NULL;
    }
    if ((*p)->x != NULL) {
        free((*p)->x);
        (*p)->x = NULL;
    }
    if ((*p)->min != NULL) {
        free((*p)->min);
        (*p)->min = NULL;
    }
    if ((*p)->max != NULL) {
        free((*p)->max);
        (*p)->max = NULL;
    }
    free(*p);
    *p = NULL;
}

void bound(node_t *p, set_t **h, double *zp, double *x) {
    if (p == NULL || h == NULL) {
        return;
    }
    if(p->z > *zp) {
        *zp = p->z;
        memcpy(x, p->x, (p->n) * sizeof(double));

        if ((*h) == NULL) {
            return;
        }
        set_t *q, *prev, *next;
        q = *h;
        while(q->data->z < p->z) {
            q = q->next;
            if(q == NULL) {
                return;
            }
            if(q->data == NULL) {
                return;
            }
        }

        prev = q;

        q = q->next;
        while(q != NULL) {
            next = q->next;
            if (q->data->z < p->z) {
                prev->next = q->next;
                free_node(&q->data);
                free(q);
            } else {
                prev = q;
            }
            q = next;
        }
    }
}

int branch(node_t *q, double z) {
    double min,max;
    int h;
    if(q->z < z) {
        return 0;
    }

    for(h=0; h < q->n; h++) {
        if(!is_integer(&q->x[h])) {
            if(q->min[h] == -INFINITY) {
                min = 0;
            } else {
                min = q->min[h];
            }
            max = q->max[h];
            if ((floorf(q->x[h]) < min) || (ceilf(q->x[h] > max))) {
                continue;
            }
            q->h = h;
            q->xh = q->x[h];
            
            // delete each of a,b,c,x
            if (q->a != NULL) { // delete matrix
                for (int i = 0; i < q->m + 1; i++) {
                    free(q->a[i]);
                    q->a[i] = NULL;
                }
                free(q->a);
                q->a = NULL;
            }
            if (q->b != NULL) { // etc..
                free(q->b);
                q->b = NULL;
            }
            if (q->c != NULL) {
                free(q->c);
                q->c = NULL;
            }
            if (q->x != NULL) {
                free(q->x);
                q->x = NULL;
            }
            return 1;
        }
    }
    return 0;
}



void succ(node_t *p, set_t **h, int m, int n, double **a, double *b, double *c, int k, double ak, double bk, double *zp, double *x) {
    node_t *q = extend(p, m, n, a, b, c, k, ak, bk);

    if(q==NULL) {
        return;
    }

    q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);
    if(isfinite(q->z)) {
        if (integer(q)) {
            bound(q, h, zp, x);
            free_node(&q);
        } else if(branch(q, *zp)) {
            set_t *node = malloc(sizeof(set_t));
            node->data = q;
            node->next = NULL;

            //Insert new node
            if(h != NULL) {
                node->next = *h;
                *h = node;
            } else {
                *h = node;
            }
            return;
        }
    }
    free_node(&q);
}

void free_set(set_t **h) {
    set_t *current = *h;
    set_t *next;

    while(current != NULL) {
        next = current->next;
        free_node(&current->data);
        free(current);
        current = next;
    }
    *h = NULL;
}

double intopt(int m, int n, double **a, double *b, double *c, double *x) {
    node_t* p = initial_node(m, n, a, b, c);
    set_t* h = calloc(m, sizeof(set_t));

    h->data = p;
    h->next = NULL;

    double z = -INFINITY;

    p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);

    if ((integer(p)) || !isfinite(p->z)) {
        z = p->z;
        if (integer(p)) {
            memcpy(x, p->x, (p->n)*sizeof(double));
        }
        free_node(&p);
        free(h);
        return z;
    }
    branch(p,z);
    while(h != NULL) {
        set_t *pop = h;
        h = pop->next;
        node_t *p = pop->data;
        free(pop);
        succ(p, &h, m, n, a, b, c, p->h, 1, floorf(p->xh), &z, x);
        succ(p, &h, m, n, a, b, c, p->h, -1, -ceilf(p->xh), &z, x);
        free_node(&p);
    }

    free_set(&h);
    h = NULL;

    if (z == -INFINITY) {
        return NAN;
    } else {
        return z;
    }

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

void print_matrix(double** a, int m, int n) {
    // Print matrix A
    printf("Matrix A (m x (n + 2)):\n");
    for (int i = 0; i < m; i++) {
        printf("Row %d: ", i);
        for (int j = 0; j < n + 2; j++) {
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
}

void print_simplex(simplex_t *s) {
    if (s == NULL) {
        printf("Simplex structure is NULL.\n");
        return;
    }

    printf("Simplex structure:\n");
    printf("Number of coefficients (m): %d\n", s->m);
    printf("Number of decision variables (n): %d\n", s->n);

    // Print var array
    //printf("Variable indices (var):\n");
    //for (int i = 0; i < s->n + s->m + 1; i++) {
    //    printf("var[%d]: %d\n", i, s->var[i]);
    //}

    // Print matrix A
    printf("Matrix A (m x (n + 2)):\n");
    for (int i = 0; i < s->m; i++) {
        printf("Row %d: ", i);
        for (int j = 0; j < s->n + 2; j++) {
            printf("%lf ", s->a[i][j]);
        }
        printf("\n");
    }

    // Print vector b
    printf("Vector b:\n");
    for (int i = 0; i < s->m; i++) {
        printf("b[%d]: %lf\n", i, s->b[i]);
    }

    // Print vector x
    //printf("Vector x:\n");
    //for (int i = 0; i < s->n + 1; i++) {
    //    printf("x[%d]: %lf\n", i, s->x[i]);
    //}

    // Print vector c
    printf("Vector c:\n");
    for (int i = 0; i < s->n; i++) {
        printf("c[%d]: %lf\n", i, s->c[i]);
    }

    // Print objective function value
    printf("Objective function value (y): %lf\n", s->y);
}


int main() {
    int m, n;
    double coeff, matrix_value, constants;
    double result;
    double y = 0.0;
    //printf("Enter the size of matrix and  vectors: \n");
    scanf("%d %d", &m, &n);

    double** a = make_matrix(m,n);

    double* b = make_constants(m);

    double* c = make_coeff_vector(n);

    double* x = make_dec_variables(n);

    for(int j = 0; j < n; j+=1) {
        //printf("Enter coeff %d:\n", j);
        scanf("%lf", &coeff);
        c[j] = coeff;
    }

    for(int i = 0; i < m; i+=1) {
        for(int j = 0; j< n; j+=1) {
            //printf("Enter matrix_value(%d,%d):\n", i, j);
            scanf("%lf", &matrix_value);
            a[i][j] = matrix_value;
        }
    }

    for(int i = 0; i < m; i+=1) {
        //printf("Enter dec variable %d:\n", i);
        scanf("%lf", &constants);
        b[i] = constants;
    }

   // Print the matrix and vectors
    //print_coeff_equation(c, n);
    //print_inequalities(a, b, m, n);

    result = intopt(m, n, a, b, c, x);
    
    printf("Output: %10.3lf\n", result);
    // Free allocated memory
    for (int i = 0; i < m; i++) {
        free(a[i]);
    }
    free(a);
    free(c);
    free(b);
    free(x);

    return 0;
}
