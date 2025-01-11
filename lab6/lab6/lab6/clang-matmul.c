/* Anything may be changed, except matrix sizes and float type. 
 *
 * The output when defining CHECK should be the same with and without tuning.
 * 
 */

#include <stdio.h>
#include <string.h>

#define N (512)

#ifdef CHECK
#define TYPE	int
#else
#define TYPE	float
#endif

TYPE a[N][N];
TYPE b[N][N];
TYPE c[N][N];

void matmul() {
  size_t i, j, k;
  
  for (i = 0; i < N; i++) {
    for (k = 0; k < N; k++) {
      TYPE bik = b[i][k];  // Pull this out to help vectorization
      for (j = 0; j < N; j++) {
        a[i][j] += bik * c[k][j];
      }
    }
  }
}

void init()
{
#ifdef CHECK
	size_t	i, j;

	for (i = 0; i < N; i += 1) {
		for (j = 0; j < N; j += 1) {
			b[i][j] = 12 + i * j * 13;
			c[i][j] = -13 + i + j * 21;
		}
	}
#endif
}

void output()
{
#ifdef CHECK
	size_t	i, j;

	for (i = 0; i < N; i += 1)
		for (j = 0; j < N; j += 1)
			printf("a[%3zu][%3zu] = %d\n", i, j, a[i][j]);
#endif
}

int main()
{
	init();
	matmul();
	output();

	return 0;
}
