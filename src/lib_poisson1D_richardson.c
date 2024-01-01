/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <math.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

void eig_poisson1D(double *eigval, int *la)
{
  double B = 2;

  for (int i = 1; i < *la + 1; i++)
  {
    eigval[i - 1] = B + 2 * (-1) * cos(i * M_PI / (*la + 1));
  }
  return;
}

double eigmax_poisson1D(int *la)
{
  double A = -1;
  double B = 2;
  double C = -1;
  int K = *la;

  return B + 2 * (-1) * cos((*la) * M_PI / (*la + 1));
}

double eigmin_poisson1D(int *la)
{

  double A = -1;
  double B = 2;
  double C = -1;
  int K = 1;

  return B + 2 * (-1) * cos(K * M_PI / (*la + 1));
}

double richardson_alpha_opt(int *la)
{
  double eig_min = eigmin_poisson1D(la);
  double eig_max = eigmax_poisson1D(la);
  return 2 / (eig_min + eig_max);
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  double *Y = malloc(*la * sizeof(double));
  double norme_b = cblas_dnrm2(*la, RHS, 1);

  for ((*nbite) = 0; (*nbite) < *maxit; (*nbite)++)
  {
    // Copie de RHS => Y
    cblas_dcopy(*la, RHS, 1, Y, 1);

    // Y = b - Ax
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, Y, 1);

    // x = x + alpha*(b-Ax)
    // x = x + alpha * Y
    cblas_daxpy(*la, *alpha_rich, Y, 1, X, 1);

    // Calcul du residu
    resvec[*nbite] = cblas_dnrm2(*la, Y, 1) / norme_b;

    if (resvec[*nbite] <= *tol)
      break;
  }
  free(Y);
  return;
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv)
{
  int iterator = 0;
  for (int i = 0; i < (*la); i++)
  {
    MB[iterator + *kv] = AB[iterator + *kv];
    iterator += *lab;
  }
  return;
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv)
{
  int iterator = 0;
  for (int i = 0; i < (*la); i++)
  {
    MB[iterator + *kl] = AB[iterator + *kl];
    MB[iterator + 1 + *kl] = AB[iterator + 1 + *kl];
    iterator += *lab;
  }
  return;
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite)
{

  double *Y = malloc(*la * sizeof(double));
  double norme_b = cblas_dnrm2(*la, RHS, 1);

  int *ipiv = malloc(*la * sizeof(int));
  int info = 0;
  int NRHS = 1;
  int ku_moins1 = *ku - 1;

  // Factorisation LU
  dgbtrf_(la, la, kl, &ku_moins1, MB, lab, ipiv, &info);

  for ((*nbite) = 0; (*nbite) < *maxit; (*nbite)++)
  {
    // Copie de RHS => Y
    cblas_dcopy(*la, RHS, 1, Y, 1);

    // b = b - Ax
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, Y, 1);

    // Calcul du residu
    resvec[*nbite] = cblas_dnrm2(*la, Y, 1) / norme_b;

    // b = b/M | b = (b - Ax)/M
    dgbtrs_("N", la, kl, &ku_moins1, &NRHS, MB, lab, ipiv, Y, la, &info);

    // x = x + b | x = x + (b - Ax)/M
    cblas_daxpy(*la, 1, Y, 1, X, 1);
    if (resvec[*nbite] <= *tol)
      break;
  }
  free(Y);
  free(ipiv);

  return;
}
