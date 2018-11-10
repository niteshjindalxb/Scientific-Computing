#include <bits/stdc++.h>
#include "main.h"
#include "parabolic_forward_diff.h"
using namespace std;

// -----------------------------------------------------------
double forward_stencil (double x, double t, std::vector<double> prev_soln)
{
    double x_0 = x_range_from();
    double delta_x = x_step_size();
    double delta_t = t_step_size();
    double c = parabolic_c();
    double lambda = (c*c*delta_t)/(delta_x*delta_x);

    // Calculate indices
    int i = (x-x_0)/delta_x;
    double soln = lambda*prev_soln[i+1] + (1.0-2.0*lambda)*prev_soln[i] + lambda*prev_soln[i-1];
    return soln;
}
// -----------------------------------------------------------
std::vector<double> forward_diff_method ()
{
    double x_0 = x_range_from();
    double x_1 = x_range_to();
    double t_0 = t_range_from();
    double t_1 = evaluate_at();
    double delta_x = x_step_size();
    double delta_t = t_step_size();

    std::vector<double> prev_soln = initial_value ();
    int t_size = (t_1 - t_0)/delta_t + 1.5;
    int x_size = (x_1 - x_0)/delta_x + 1.5;

    for (int t_h = 1; t_h < t_size; t_h++)
    {
        double t = t_0 + delta_t * t_h;
        std::vector<double> cur_soln;
        // Find soln at t
        for (int x_h = 0; x_h < x_size; x_h++)
        {
            double x = x_0 + delta_x * x_h;
            // if x == x_0 or x == x_1
            // then use boundary_condition
            if (x_h == 0 || x_h == x_size - 1)
                cur_soln.push_back (boundary_cond());
            else
            // Calculate value at current point (x, t) using stencil
                cur_soln.push_back (forward_stencil(x, t, prev_soln));
        }
        prev_soln = cur_soln;
    }
    return prev_soln;
}
