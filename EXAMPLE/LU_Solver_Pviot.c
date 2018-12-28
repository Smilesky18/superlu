# include "slu_ddefs.h"
# include <stdio.h>

void LU(double A[][800], double x[800], int n)
{
  double L[n][n], U[n][n];
  double sum_U = 0.0, sum_L = 0.0, sum_y = 0.0, sum_x = 0.0;
  double y[n];
  double check_sum[n], sum, divide;
  int i, j, r, k, check;
  int record_order[n];
  double temp = 0.0;
  double max = 0.0;
  int record_order_temp;
  FILE *fp_L = fopen("L_pviot.txt", "w");
  FILE *fp_U = fopen("U_pviot.txt", "w");
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
	L[i][j] = 0.0;
	U[i][j] = 0.0;
    }
  }
  for ( i = 0; i < n; i++ )
  {
    record_order[i] = i;
  }
  for ( i = 0; i < n; i++ )
  {
    if ( temp < A[i][0] )
    {
      record_order_temp = i;
      temp = A[i][0];
    }
  }
  if ( record_order_temp == record_order[0] )
  {
    printf("do not need to change the pviot!");
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      temp = A[0][i];
      A[0][i] = A[record_order_temp][i];
      A[record_order_temp][i] = temp;
      //printf("the %dth time to change the pviot!", record_order_temp);
    }
    temp = 0.0;
  }
  // 将调换的次序记录在record_order中，然后调换矩阵中的数值
  record_order[0] = record_order_temp;
  
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
  // 第二列以及以下的元素置为零
  for ( i = 1; i < n; i++ )
  {
    if ( A[i][0] == 0 )
      continue;
    else
    {
      divide = A[i][0]/A[0][0];
      for ( j = 0; j < n; j++ )
      {
	A[i][j] = A[i][j] - divide * A[0][j];
      }
    } 
  }
    
  // 从第1步开始，进行选主元的LU分解法
  for ( r = 1; r < n; r++ )
  {
    for ( k = 0; k < r-1; k++ )
    {
      sum += L[r][k]*U[k][r];
    }
    for ( i = r; i < n; i++ )
    {
      
      check_sum[r] = A[i][r] - sum;
      //sum = 0.0;
      if ( check_sum[r] > max )
      {
	max = check_sum[r];
	record_order_temp = i;
      }
    }
   // sum = 0.0;
   // max = 0.0;
    if ( record_order_temp == record_order[r] )
    {
      printf (" do not need to change ");
      sum = 0.0;
    }
    else
    {
      for ( i = 0; i < n; i++ )
      {
	temp = A[r][i];
	//printf("record_order[%d] is: %d", r, record_order[r]);
	A[r][i] = A[record_order_temp][i];
	A[record_order_temp][i] = temp;
	record_order[r] = record_order_temp;
	//printf("the %dth time to change the pviot!", record_order_temp);
      }
    }
    // 计算U的第r行元素和L的第r列元素
    U[r][r] = A[r][r] - sum;
    sum = 0.0;
    max = 0.0;
    for ( i = r+1; i < n; i++ )
    {
      L[i][r] = A[i][r] / U[r][r];
      for ( k = 0; k < r-1; k++ )
      {
	sum_U += L[r][k]*U[k][i];
      }
      U[r][i] = A[r][i] - sum_U;
      sum_U = 0.0;
    }
    // 将r行以下的元素置为0
    for ( i = r+1; i < n; i++ )
    {
      if ( A[i][r] == 0 )
	continue;
      else
      {
	divide = A[i][r] / A[r][r];
	for ( j = r; j < n; j++ )
	{
	  A[i][j] = A[i][j] - divide * A[r][j];
	}
      }
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
    FILE *fp_A_ORI = fopen("A_ORI.txt", "w");
    FILE *fp_A_LU = fopen("A_LU.txt", "w");
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
       for ( j = 0; j < n; j++ )
       {
	  fprintf(fp_A_ORI, "%lf	",A[i][j]);
       }
       fprintf(fp_A_ORI, "\n");
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
    /*for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
	printf("A[%d][%d]=%f ", i, j, A[i][j]);
      }
    }*/
    for ( i = 0; i < n; i++ )
    {
       for ( j = 0; j < n; j++ )
       {
	  fprintf(fp_A_LU, "%lf	",A[i][j]);
       }
       fprintf(fp_A_LU, "\n");
    }
     
    /*for ( i = 0; i < n; i++ )
    {
      printf("在main函数中的解元素为: %lf\n", x[i]);
    }*/	
    return 0;
}

