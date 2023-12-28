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
  return 0;
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
