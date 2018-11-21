/* Solve the following parabolic partial differential equation */
/* Elliptic partial differential equation */

#include <bits/stdc++.h>
#include "main.h"
#include "gauss_seidel.h"
#include "display.h"
using namespace std;

// -----------------------------------------------------------
double orig_sol (double x, double y)
{
    return 100.0 * sinh(M_PI*y/10.0) * sin(M_PI*x/10.0) / sinh(M_PI*15/10.0);
}
/* Boundary values */
double boundary_cond (double x, double y)
{
    if (abs(x-x_range_from) < epsilon || abs(x-x_range_to) < epsilon || abs(y-y_range_from) < epsilon)
        return 0.0;
    else /* y == 15 */
        return 100.0 * sin(M_PI*x/10.0);
}
std::vector<double> calc_exact_soln ()
{
    vector<double> exact_soln;
    // Calculate the mesh size
    double h = (x_range_to - x_range_from)/x_grid_size;
    double k = (y_range_to - y_range_from)/y_grid_size;

    for (int j = y_grid_size-1; j >= 0; j--)
        for (int i = 1; i < x_grid_size; i++)
        {
            double x = x_range_from + i*h;
            double y = y_range_from + j*k;
            exact_soln.push_back (orig_sol(x, y));
        }
    return exact_soln;
}
void init_boundary (std::vector<std::vector<double> > &grid_value)
{
    // Calculate the mesh size
    double h = (x_range_to - x_range_from)/x_grid_size;
    double k = (y_range_to - y_range_from)/y_grid_size;

    // Lower horizontal boundary
    for (int i = 0; i < x_grid_size+1; i++)
    {
        int j = y_grid_size;
        double x = x_range_from + i*h;
        double y = y_range_to - j*k;
        grid_value[j][i] = boundary_cond(x, y);
    }
    // Upper horizontal boundary
    for (int i = 0; i < x_grid_size+1; i++)
    {
        int j = 0;
        double x = x_range_from + i*h;
        double y = y_range_to - j*k;
        grid_value[j][i] = boundary_cond(x, y);
    }
    // Left vertical boundary
    for (int j = 0; j < y_grid_size+1; j++)
    {
        int i = 0;
        double x = x_range_from + i*h;
        double y = y_range_to - j*k;
        grid_value[j][i] = boundary_cond(x, y);
    }
    // Right vertical boundary
    for (int j = 0; j < y_grid_size+1; j++)
    {
        int i = x_grid_size;
        double x = x_range_from + i*h;
        double y = y_range_to - j*k;
        grid_value[j][i] = boundary_cond(x, y);
    }

}
vector<double> e_five_pt_solver (std::vector<std::vector<double> > &grid_value)
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

    for (int i = 0; i < num_x; i++)
    {
        for (int j = 0; j < num_y; j++)
        {
            int k = i + j*num_x;
            // Input into row k
            coef_A[k][i + j*num_x] = b;

            if (i+2 < num_x+1)
                coef_A[k][i + 1 + j*num_x] = -1.0;
            else
                coef_b[k] += grid_value[i+2][j+1];

            if (i > 0)
                coef_A[k][i - 1 + j*num_x] = -1.0;
            else
                coef_b[k] += grid_value[i][j+1];

            if (j+2 < num_y+1)
                coef_A[k][i + (j+1)*num_x] = -a;
            else
                coef_b[k] += grid_value[i+1][j+2] * a;

            if (j > 0)
                coef_A[k][i + (j-1)*num_x] = -a;
            else
                coef_b[k] += grid_value[i+1][j] * a;
        }
    }
    for (int j=0; j<num_y+2; j++)
    {
        for (int i=0; i<num_x+2; i++)
            cout << grid_value[j][i] << " ";
        cout << endl;
    }
    display(coef_b);
    vector<double> soln = gauss_seidel(coef_A, coef_b);
    return soln;
}
// -----------------------------------------------------------
int main()
{
    std::vector<std::vector<double> > grid_value (y_grid_size+1, std::vector<double>(x_grid_size+1, 0.0));
    vector<double> soln = e_five_pt_solver (grid_value);
    vector<double> exact_soln = calc_exact_soln ();

    display_2_vectors (soln, exact_soln);
    return 0;
}
