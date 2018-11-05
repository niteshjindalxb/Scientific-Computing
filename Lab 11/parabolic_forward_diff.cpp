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

    for (double t = t_0 + delta_t; t < t_1 + epsilon; t += delta_t)
    {
        std::vector<double> cur_soln;
        // Find soln at t
        for (double x = x_0; x < x_1 + epsilon; x += delta_x)
        {
            // if x == x_0 or x == x_1
            // then use boundary_condition
            if (x < x_0 + epsilon || (x < x_1 + epsilon && x > x_1 - epsilon))
                cur_soln.push_back (boundary_cond());
            else
            // Calculate value at current point (x, t) using stencil
                cur_soln.push_back (forward_stencil(x, t, prev_soln));
        }
        prev_soln = cur_soln;
    }
    return prev_soln;
}
