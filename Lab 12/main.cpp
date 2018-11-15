/* Solve the following parabolic partial differential equation */
/* Elliptic partial differential equation */

#include <bits/stdc++.h>
#include "main.h"
// #include "gauss_seidel.h"
using namespace std;

// -----------------------------------------------------------
double orig_sol (double x, double y)
{
    return 100.0 * sinh(M_PI*y/10.0) * sin(M_PI*x/10.0) / sinh(M_PI*15/10.0);
}
/* Boundary values */
double boundary_cond (double x, double y)
{
    if (abs(x) < epsilon || abs(x-x_range_from) < epsilon || abs(y) < epsilon)
        return 0.0;
    else /* y == 15 */
        return 100.0 * sin(M_PI*x/10.0);
}
// std::vector<double> calc_exact_soln ()
// {
//     double x_0 = x_range_from;
//     double x_1 = x_range_to;
//     double delta_x = x_step_size;
//     double t = evaluate_at;
//
//     std::vector<double> exact_soln;
//
//     for (double x = 0; x < x_1 + epsilon; x += delta_x)
//         exact_soln.push_back (orig_sol(x, t));
//
//     return exact_soln;
// }
// std::vector<double> calc_error (vector<double> &approx_soln, vector<double> &exact_soln)
// {
//     vector<double> list_error;
//     for (int i = 0; i < exact_soln.size(); i++)
//         list_error.push_back (abs(approx_soln[i] - exact_soln[i]));
//     return list_error;
// }
void init_boundary (std::vector<std::vector<double> > &grid_value)
{
    // Calculate the mesh size
    double h = (x_range_to - x_range_from)/x_grid_size;
    double k = (y_range_to - y_range_from)/y_grid_size;

    for (int i = 0; i < x_grid_size+1; i++)
    {
        for (size_t j = 0; j < y_grid_size+1; j++)
        {
            if (j != 0 && j != y_grid_size)
                continue;
            double x = x_range_from + i*h;
            double y = y_range_from + j*k;
            std::cout << "i = " << i << '\n';
            std::cout << "j = " << j << '\n';
            grid_value[i][j] = boundary_cond(x, y);
        }
    }
}
void e_five_pt_solver (std::vector<std::vector<double> > &grid_value)
{
    init_boundary (grid_value);
    int num_x = grid_value.size()-2;
    int num_y = grid_value[0].size()-2;
    // Calculate the mesh size
    double h = (x_range_to - x_range_from)/x_grid_size;
    double k = (y_range_to - y_range_from)/y_grid_size;

    // Evaluate the coefficient matrix
    int nVar = num_x * num_y;
    std::vector<std::vector<double> > coef_A(nVar, std::vector<double>(nVar, 0.0));
    std::vector<double> coef_b(nVar, 0.0);

    double a = (h/k)*(h/k);
    double b = 2.0 * (1.0 + a);

    for (size_t i = 0; i < num_x; i++)
    {
        for (int j = 0; j < num_y; j++)
        {
            int k = i + j*num_x;
            // Input into row k
            coef_A[k][i + j*num_x] = b;

            if (i < num_x-1)
                coef_A[k][i + 1 + j*num_x] = -1.0;
            else
                coef_b[k] += grid_value[i+2][j+1];

            if (i > 0)
                coef_A[k][i - 1 + j*num_x] = -1.0;
            else
                coef_b[k] += grid_value[i][j+1];

            if (j < num_y-1)
                coef_A[k][i + (j+1)*num_x] = -a;
            else
                coef_b[k] += grid_value[i+1][j+2] * a;

            if (j > 0)
                coef_A[k][i + (j-1)*num_x] = -a;
            else
                coef_b[k] += grid_value[i+1][j] * a;
        }
    }
    // for (size_t i = 0; i < coef_A.size(); i++) {
    //     for (size_t j = 0; j < coef_A[0].size(); j++) {
    //         // std::cout << "[" << i << "][" << j << "]" << coef_A[i][j] << ' ';
    //         std::cout << fixed << coef_A[i][j] << ' ';
    //     }
    //     std::cout << '\n';
    // }
    for (size_t i = 0; i < coef_b.size(); i++) {
        std::cout << coef_b[i] << '\n';
    }
}
// -----------------------------------------------------------
int main()
{
    std::vector<std::vector<double> > grid_value (x_grid_size+1, std::vector<double>(y_grid_size+1, 0.0));
    e_five_pt_solver (grid_value);
    // for (size_t i = 0; i < grid_value.size(); i++) {
    //     for (size_t j = 0; j < grid_value[0].size(); j++) {
    //         std::cout << "[" << i << "][" << j << "]" << grid_value[i][j] << '\n';
    //     }
    // }
    return 0;
}
