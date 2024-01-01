/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
void set_GB_operator_colMajor_poisson1D(double *AB, int *lab, int *la, int *kv)
{
  int iterator = 0;
  for (int i = 0; i < (*la); i++)
  {
    if (i == 0)
    {
      AB[iterator + 1 + *kv] = 2.0;
      AB[iterator + 2 + *kv] = -1.0;
    }

    else if (i == (*la) - 1)
    {
      AB[iterator + *kv] = -1.0;
      AB[iterator + 1 + *kv] = 2.0;
    }

    else
    {
      AB[iterator + *kv] = -1.0;
      AB[iterator + 1 + *kv] = 2.0;
      AB[iterator + 2 + *kv] = -1.0;
    }

    iterator += *lab;
  }

  return;
}

void set_GB_operator_colMajor_poisson1D_Id(double *AB, int *lab, int *la, int *kv)
{

  int iterator = 0;
  for (int i = 0; i < (*la); i++)
  {
    AB[iterator + 1 + *kv] = 1.0;
    iterator += *lab;
  }
  return;
}

void set_dense_RHS_DBC_1D(double *RHS, int *la, double *BC0, double *BC1)
{
  RHS[0] = *BC0;
  RHS[*la - 1] = *BC1;
  for (int i = 1; i < *la - 1; ++i)
  {
    RHS[i] = 0.0;
  }
  return;
}

// Solution analytique à notre probleme : T(x) = T0 + x(T1 − T0)
void set_analytical_solution_DBC_1D(double *EX_SOL, double *X, int *la, double *BC0, double *BC1)
{
  for (int i = 0; i < (*la); i++)
  {
    EX_SOL[i] = *BC0 + X[i] * (*BC1 - *BC0);
  }
  return;
}

void set_grid_points_1D(double *x, int *la)
{
  double h = 1.0 / (*la + 1);
  for (int i = 0; i < (*la); i++)
  {
    x[i] = (i + 1) * h;
  }
  return;
}

double relative_forward_error(double *x, double *y, int *la)
{
  double norm_xy = 0;
  double norm_x = 0;
  for (int i = 0; i < *la; ++i)
  {
    double ecart = x[i] - y[i];
    norm_x += x[i] + x[i];
    norm_xy += ecart * ecart;
  }
  norm_x = sqrt(norm_x);
  norm_xy = sqrt(norm_xy);
  return norm_xy / norm_x;
}

int indexABCol(int i, int j, int *lab)
{
  return 0;
}

int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info)
{
  // ipiv le vecteur qui indique les permutations
  ipiv[0] = 1;
  info = 0;
  int iterator = 0;
  for (int i = 1; i < *la; i++)
  {
    AB[iterator - 1] /= AB[iterator - 2];
    AB[iterator + *lab - 2] -= AB[iterator - 1] * AB[iterator + *lab - 3];
    ipiv[i] = i + 1;
    iterator += *lab;
  }
  return *info;
}

void set_CSR_poisson1D(int *n, double *AA, int *JA, int *IA)
{
  // taille AA = (nb element non nul * sizeof(double));
  // taille JA = (nb element non nul * sizeof(double));
  // taille IA = ( (n+1)*sizeof(double));

  int iterator = 0;

  for (int i = 0; i < *n; ++i)
  {
    if (i == 0)
    {
      IA[i] = i;
      AA[iterator] = 2;
      AA[iterator + 1] = -1;

      JA[iterator] = i;
      JA[iterator + 1] = i + 1;

      iterator += 2;
    }
    else if (i == *n - 1)
    {
      IA[i] = i;
      AA[iterator] = -1;
      AA[iterator + 1] = 2;

      JA[iterator] = i;
      JA[iterator + 1] = i + 1;

      iterator += 2;
    }
    else
    {
      IA[i] = i;
      AA[iterator] = -1;
      AA[iterator + 1] = 2;
      AA[iterator + 2] = -1;

      JA[iterator] = i - 1;
      JA[iterator + 1] = i;
      JA[iterator + 2] = i + 1;

      iterator += 3;
    }
  }
  return;
}

void set_CSC_poisson1D(int *n, double *AA, int *JA, int *IA)
{
  // AA = malloc(nb element non nul * sizeof(double));
  // JA = malloc(nb element non nul * sizeof(double));
  // IA = malloc( (n+1)*sizeof(double));

  int iterator = 0;

  for (int i = 0; i < *n; ++i)
  {
    if (i == 0)
    {
      IA[i] = iterator;
      AA[iterator] = 2;
      AA[iterator + 1] = -1;

      JA[iterator] = i;
      JA[iterator + 1] = i + 1;

      iterator += 2;
    }
    else if (i == *n - 1)
    {
      IA[i] = iterator;
      AA[iterator] = -1;
      AA[iterator + 1] = 2;

      JA[iterator] = i;
      JA[iterator + 1] = i + 1;

      iterator += 2;
    }
    else
    {
      IA[i] = iterator;
      AA[iterator] = -1;
      AA[iterator + 1] = 2;
      AA[iterator + 2] = -1;

      JA[iterator] = i - 1;
      JA[iterator + 1] = i;
      JA[iterator + 2] = i + 1;

      iterator += 3;
    }
  }
  return;
}

double *dcsrmv(int *n, double *AA, int *JA, int *IA, double *vec)
{
  double *res = malloc(*n * sizeof(double));
  for (int i = 0; i < *n; i++)
  {
    res[i] = 0.0;
    for (int j = IA[i]; j < IA[i + 1]; j++)
    {
      res[i] += AA[j] * vec[JA[j]];
    }
  }
  return res;
}


double *dcscmv(int *n, double *AA, int *JA, int *IA, double *vec)
{
  double *res = malloc(*n * sizeof(double));
  for (int i = 0; i < *n; i++)
  {
    res[i] = 0.0;
    for (int j = IA[i]; j < IA[i + 1]; j++)
    {
      res[i] += AA[j] * vec[JA[j]];
    }
  }
  return res;
}