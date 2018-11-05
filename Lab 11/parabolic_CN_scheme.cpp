#include <bits/stdc++.h>
#include "main.h"
#include "parabolic_CN_scheme.h"
#include "gauss_seidel.h"
using namespace std;

// -----------------------------------------------------------
std::vector<double> vector_stencil (std::vector<double> prev_soln, double lambda)
{
    std::vector<double> new_arr(prev_soln);
    for (size_t i = 1; i < prev_soln.size()-1; i++)
        new_arr[i] = lambda*prev_soln[i-1] + (2 - 2*lambda)*prev_soln[i] + lambda*prev_soln[i+1];

    return new_arr;
}
// -----------------------------------------------------------
std::vector<std::vector<double> > CN_coef_matrix (double lambda)
{
    double x_0 = x_range_from();
    double x_1 = x_range_to();
    double delta_x = x_step_size();
    int size = (x_1 + epsilon - x_0)/delta_x + 1;

    std::vector<std::vector<double> > matrix(size, std::vector<double> (size, 0));
    for (int i = 0; i < size; i++)
    {
        if (i == 0 || i == size-1)
            matrix[i][i] = 1;
        else
        {
            matrix[i][i-1] = -lambda;
            matrix[i][i]   = 2 + 2*lambda;
            matrix[i][i+1] = -lambda;
        }
    }
    return matrix;
}

// -----------------------------------------------------------
std::vector<double> CN_method ()
{
    double x_0 = x_range_from();
    double x_1 = x_range_to();
    double t_0 = t_range_from();
    double t_1 = evaluate_at();
    double delta_x = x_step_size();
    double delta_t = t_step_size();
    double c = parabolic_c();
    double lambda = (c*c*delta_t)/(delta_x*delta_x);

    std::vector<double> prev_soln = initial_value ();

    // Make matrix A
    std::vector<std::vector<double> > A = CN_coef_matrix(lambda);

    // Find soln at time t
    for (double t = t_0 + delta_t; t < t_1 + epsilon; t += delta_t)
    {
        prev_soln = vector_stencil (prev_soln, lambda);
        std::vector<double> cur_soln = gauss_seidel(A, prev_soln);
        if (cur_soln.size() != prev_soln.size())
            print(t);
        prev_soln = cur_soln;
    }

    return prev_soln;
}
