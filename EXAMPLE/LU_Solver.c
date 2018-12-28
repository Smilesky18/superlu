# include "slu_ddefs.h"
# include <stdio.h>

void LU(double A[][800], double x[800], int n)
{
  double L[n][n], U[n][n];
  double sum_U = 0.0, sum_L = 0.0, sum_y = 0.0, sum_x = 0.0;
  double y[n];
  int i, j, r, k;
  FILE *fp_L = fopen("L_solver.txt", "w");
  FILE *fp_U = fopen("U_solver.txt", "w");
  //U的第一行的解向量
  for ( i = 0; i < n; i++ )
  {
    U[0][i] = A[0][i];
  }
  //L的第一列的解向量
  for ( i = 1; i < n; i++ )
  {
    L[i][0] = A[i][0]/U[0][0];
  }
  //对L的对角线元素赋值为1
  for ( i = 0; i < n; i++ )
  {
    L[i][i] = 1;
  }
  //计算U的第r行，L的第r列元素（r=2,...n-1）
  for ( r = 1; r < n; r++ )
  {
    for ( i = r; i < n; i++ )
    {
      for ( k = 0; k < r; k++ )
      {
	sum_U += L[r][k]*U[k][i];
      }
      U[r][i] = A[r][i] - sum_U;
      sum_U = 0.0;
    }
    for ( i = r+1; i < n; i++ )
    {
      for ( k = 0; k < r; k++ )
      {
	sum_L += L[i][k]*U[k][r];
      }
      L[i][r] = ( A[i][r] - sum_L ) / U[r][r];
      sum_L = 0.0;
    }
  }
  printf("The element of L: ");
  for ( i = 0; i < n; i++ )
    for( j = 0; j < n; j++ )
    {
      //printf("L[%d][%d]=%f ", i, j, L[i][j]);
      fprintf(fp_L, "L[%d][%d]=%lf\n", i, j, L[i][j]);
    }
  
  printf("The element of U: ");
  for ( i = 0; i < n; i++ )
    for( j = 0; j < n; j++ )
    {
      //printf("U[%d][%d]=%f ", i, j, U[i][j]);
      fprintf(fp_U, "U[%d][%d]=%lf\n", i, j, U[i][j]);
    }
  //求解Ly=b, Ux=y的计算公式
  y[0] = A[0][n];  
  for ( i = 1; i < n; i++ )
  {
    for ( k = 0; k < i; k++ )
    {
      sum_y += L[i][k]*y[k]; 
    }
    y[i] = A[i][n] - sum_y;
    sum_y = 0.0;
  }
  x[n-1] = y[n-1] / U[n-1][n-1];
  for ( i = n-2; i >= 0; i-- )
  {
    for ( k = i + 1; k < n; k++ )
    {
      sum_x += U[i][k]*x[k];
    }
    x[i] = ( y[i] - sum_x ) / U[i][i];
    sum_x = 0.0;
  }
  for ( i = 0; i < n; i++ )
  {
    printf("解元素为: %lf\n", x[i]);
  }
}

int main( int argc, char *argv[] )
{
    FILE *fp = stdin;
    double   *a;
    int      *asub, *xa;
    int m, n, nnz;
    //const int numarray;
    int i, j;
    int NUMCol, counter_num = 0;
    double A[800][800];
    double x[800];
    dreadMM(fp, &m, &n, &nnz, &a, &asub, &xa);
    /*for ( i = 0; i < nnz; i++ )
    {
      printf("%d\n", asub[i]);
    }*/
    printf("\n%d %d %d", m, n, nnz);
    for ( j = 0; j < n; j++ )
    {
      NUMCol = xa[j+1] - xa[j];
      for ( i = 0; i < NUMCol; i++ )
      {
	A[asub[counter_num]][j] = a[counter_num];
	counter_num++;
      }
    }
    for ( i = 0; i < n; i++ )
    {
      A[i][n] = 1.0;
    }
    /*for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
	printf("A[%d][%d]=%f ", i, j, A[i][j]);
      }
    }*/
    LU( A, x, n);
   // for ( i = 0; i < n; i++ )
   // {
   //   printf("在main函数中的解元素为: %lf\n", x[i]);
   // }	
    return 0;
}

