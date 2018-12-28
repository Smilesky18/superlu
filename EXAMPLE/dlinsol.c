/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
#include "slu_ddefs.h"

/*void LU(double A[][800], double *x[800], int n)
{
  double L[n][n], U[n][n];
  double sum_U = 0.0, sum_L = 0.0, sum_y = 0.0, sum_x = 0.0;
  double y[n];
  int i, j, r, k;
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
  for ( i = 0; i < n; i++ )LU( A, x, n);
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
  /*printf("The element of L: ");
  for ( i = 0; i < n; i++ )
    for( j = 0; j < n; j++ )
    {
      printf("L[%d][%d]=%f ", i, j, L[i][j]);
    }
  
  printf("The element of U: ");
  for ( i = 0; i < n; i++ )
    for( j = 0; j < n; j++ )
    {
      printf("U[%d][%d]=%f ", i, j, U[i][j]);
    }*/
  //求解Ly=b, Ux=y的计算公式
  /*y[0] = A[0][n];  
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
  /*for ( i = 0; i < n; i++ )
  {
    printf("解元素为: %lf\n", x[i]);
  }*/
//}

int main(int argc, char *argv[])
{
    SuperMatrix A;
    NCformat *Astore;
    DNformat *Bstore;
    int *colptr_B;
    int *rowind_B;
    double *nzval_B, *nzval_A;
    double   *a;
    int      *asub, *xa;
    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
    SuperMatrix L;      /* factor L */
    SCformat *Lstore;
    SuperMatrix U;      /* factor U */
    NCformat *Ustore;
    SuperMatrix B;
    int      nrhs, ldx, info, m, n, nnz;
    double   *xact, *rhs;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    FILE      *fp = stdin;
    int i, j;
    int NUMCol, counter_num = 0;
    double A_Me[800][800];
    double x[800];
    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif

    /* Set the default input options:
	options.Fact = DOFACT;
        options.Equil = YES;
    	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES;
     */
    set_default_options(&options);

#if 0
    /* Read the matrix in Harwell-Boeing format. */
    dreadhb(fp, &m, &n, &nnz, &a, &asub, &xa);
#else
    /* Read the matrix in Matrix Market format. */
    dreadMM(fp, &m, &n, &nnz, &a, &asub, &xa);
#endif

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    Astore = A.Store;
    nzval_A = Astore->nzval;
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
    
    nrhs   = 1;
    if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
    for ( i = 0; i < m; i++ )
    {
      rhs[i] = 1.0;
    }
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
    /*for ( i = 0; i < m; i++ )
      printf(" the element in rhs is: %f\n", rhs[i]);*/
    Bstore = (DNformat *)B.Store;
    nzval_B = (double *)Bstore->nzval;
    /*for (i=0; i < Bstore->lda; i++)
      printf("B的第%d个元素是: %lf\n",i, *(nzval_B+i));*/
    //xact = doubleMalloc(n * nrhs);
    //ldx = n;
    //dGenXtrue(n, nrhs, xact, ldx);
    //dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);

    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);
    Bstore = (DNformat *)B.Store;
    nzval_B = (double *)Bstore->nzval;
    /*for (i=0; i < Bstore->lda; i++)
      printf("before dgssv B的第%d个元素是: %lf\n",i, *(nzval_B+i));*/
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    Bstore = (DNformat *)B.Store;
    nzval_B = (double *)Bstore->nzval;
    //printf("after dgssv the lda of B is: %d\n", Bstore->lda);
    //for (i=0; i < Bstore->lda; i++)
    //  printf("after dgssv B的第%d个元素是: %lf\n",i, *(nzval_B+i));
    double *sol = (double*) ((DNformat*) B.Store)->nzval;
    for (i=0; i < Bstore->lda; i++)
      printf("解元素为: %lf\n", sol[i]);
  
    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);

    SUPERLU_FREE (rhs);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit main()");
#endif
     for ( j = 0; j < n; j++ )
    {
      NUMCol = xa[j+1] - xa[j];
      for ( i = 0; i < NUMCol; i++ )
      {
	A_Me[asub[counter_num]][j] = a[counter_num];
	counter_num++;
      }
    }
    for ( i = 0; i < n; i++ )
    {
      A_Me[i][n] = 1.0;
    }
    //LU( A, x, n);
    
}

