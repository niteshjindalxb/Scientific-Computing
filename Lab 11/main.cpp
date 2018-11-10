/* Solve the following parabolic partial differential equation */
/* du/dt - (4/pi^2)d2u/dx2 = 0 */
#include <bits/stdc++.h>
#include "main.h"
#include "parabolic_forward_diff.h"
#include "parabolic_backward_diff.h"
#include "parabolic_CN_scheme.h"
#include "gauss_seidel.h"
using namespace std;

// -----------------------------------------------------------
double orig_sol (double x, double t)
{
    return exp(-t)*sin(M_PI*x/2.0) + exp(-t/4.0)*sin(M_PI*x/4.0);
}
double x_step_size ()
{
    return 0.2;
}
double t_step_size ()
{
    return 0.04;
}
double x_range_from()
{
    return 0.0;
}
double x_range_to()
{
    return 4.0;
}
double t_range_from()
{
    return 0.0;
}
double evaluate_at()
{
    return /* at t = */ 0.4;
}
double parabolic_c()
{
    return (double) 2.0/M_PI;
}
/* Boundary values */
double boundary_cond ()
{
    return 0.0;
}
/* Initial condition */
double initial_cond (double x)
{
    return sin(M_PI*x/4.0) * (1.0 + 2.0*cos(M_PI*x/4.0));
}
/* Initial value */
std::vector<double> initial_value ()
{
    double x_0 = x_range_from();
    double x_1 = x_range_to();
    double delta_x = x_step_size();

    std::vector<double> init;
    for (double x = x_0; x < x_1 + epsilon; x += delta_x)
        init.push_back (initial_cond(x));

    return init;
}
std::vector<double> mesh_points ()
{
    std::vector<double> x_points;
    double x_0 = x_range_from();
    double x_1 = x_range_to();
    double delta_x = x_step_size();

    for (double x = x_0; x < x_1 + epsilon; x += delta_x)
        x_points.push_back (x);
    return x_points;
}
void display (vector <double> &arr)
{
    cout << fixed << setprecision(10);
    for (int i = 0; i < arr.size(); i++)
        cout << arr[i] << "\n";
    cout << endl;
}
void display_2_vectors (vector <double> &arr1, vector <double> &arr2)
{
    cout << fixed << setprecision(10);
    for (int i = 0; i < arr1.size(); i++)
        cout << arr1[i] << "\t" << arr2[i] << endl;
}
void display_5_vectors (vector <double> &arr1, vector <double> &arr2, vector <double> &arr3, vector <double> &arr4, vector <double> &arr5)
{
    cout << fixed << setprecision(10);
    for (int i = 0; i < arr1.size(); i++)
        cout << arr1[i] << "\t" << arr2[i] << "\t" << arr3[i] << "\t" << arr4[i] << "\t" << arr5[i] << endl;
}
std::vector<double> calc_exact_soln ()
{
    double x_0 = x_range_from();
    double x_1 = x_range_to();
    double delta_x = x_step_size();
    double t = evaluate_at();

    std::vector<double> exact_soln;

    for (double x = 0; x < x_1 + epsilon; x += delta_x)
        exact_soln.push_back (orig_sol(x, t));

    return exact_soln;
}
std::vector<double> calc_error (vector<double> &approx_soln, vector<double> &exact_soln)
{
    vector<double> list_error;
    for (int i = 0; i < exact_soln.size(); i++)
        list_error.push_back (abs(approx_soln[i] - exact_soln[i]));
    return list_error;
}
// -----------------------------------------------------------
int main()
{
    // Mesh x_points
    std::vector<double> x_points = mesh_points();

    // Forward difference method
    std::vector<double> forward_soln = forward_diff_method ();

    // Backward difference method
    std::vector<double> backward_soln = backward_diff_method ();

    // Crank-Nicolson method
    std::vector<double> CN_soln = CN_method ();

    // Calculate exact_soln at t = evaluate_at()
    std::vector<double> exact_soln = calc_exact_soln();

    print("Forward_soln\tExact_soln");
    display_2_vectors(forward_soln, exact_soln);

    print("--------------------------------");
    print("Backward_soln\tExact_soln");
    display_2_vectors(backward_soln, exact_soln);

    print("--------------------------------");
    print("CN_soln\tExact_soln");
    display_2_vectors(CN_soln, exact_soln);

    // Calculation of error
    std::vector<double> forward_error = calc_error (forward_soln, exact_soln);
    std::vector<double> backward_error = calc_error (backward_soln, exact_soln);
    std::vector<double> CN_error = calc_error (CN_soln, exact_soln);

    print("Forward_error");
    display(forward_error);

    print("--------------------------------");
    print("Backward_error");
    display(backward_error);

    print("--------------------------------");
    print("CN_error");
    display(CN_error);

    // display_5_vectors (x_points, forward_soln, backward_soln, CN_soln, exact_soln);
    return 0;
}
// -----------------------------------------------------------
/*
Use this to plot (using gnuplot):

plot "out.txt" using 1:2 title "Forward difference method" with lines linestyle 1, "out.txt" using 1:3 title "Backward difference method" with lines linestyle 2, "out.txt" using 1:4 title "Crank-Nicolson scheme" with lines linestyle 3, "out.txt" using 1:5 title "Exact solution" with lines linestyle 4

*/
