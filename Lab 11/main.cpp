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
    // return (x_step_size() * x_step_size() * STABILITY * M_PI * M_PI)/4.0;
    return 0.08;
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
    return /* at t = */ 4.0;
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
void display_5_vectors (vector <double> &arr1, vector <double> &arr2, vector <double> &arr3, vector <double> &arr4, vector <double> &arr5, ofstream &file_out)
{
    file_out << fixed << setprecision(10);
    for (int i = 0; i < arr1.size(); i++)
        file_out << arr1[i] << "\t" << arr2[i] << "\t" << arr3[i] << "\t" << arr4[i] << "\t" << arr5[i] << endl;
}
void display_4_vectors (vector <double> &arr1, vector <double> &arr2, vector <double> &arr3, vector <double> &arr4, ofstream &file_out)
{
    file_out << fixed << setprecision(10);
    for (int i = 0; i < arr1.size(); i++)
        file_out << arr1[i] << "\t" << arr2[i] << "\t" << arr3[i] << "\t" << arr4[i] << endl;
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
vector<double> calc_error (vector<double> &approx_soln, vector<double> &exact_soln)
{
    vector<double> list_error;
    for (int i = 0; i < exact_soln.size(); i++)
        list_error.push_back (abs(approx_soln[i] - exact_soln[i]));
    return list_error;
}
// -----------------------------------------------------------
int main()
{
    print (t_step_size());
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

    // print("Forward_soln\tExact_soln");
    // display_2_vectors(forward_soln, exact_soln);

    // print("--------------------------------");
    // print("Backward_soln\tExact_soln");
    // display_2_vectors(backward_soln, exact_soln);

    // print("--------------------------------");
    // print("CN_soln\tExact_soln");
    // display_2_vectors(CN_soln, exact_soln);

    ofstream file_out ("output_soln.txt", std::ios::out);
    display_5_vectors (x_points, forward_soln, backward_soln, CN_soln, exact_soln, file_out);
    file_out.close();

    vector<double> error_forward = calc_error (forward_soln, exact_soln);
    vector<double> error_backward = calc_error (backward_soln, exact_soln);
    vector<double> error_CN = calc_error (CN_soln, exact_soln);
    print ("\n");

    ofstream file_out2 ("output_error.txt", std::ios::out);
    display_4_vectors (x_points, error_forward, error_backward, error_CN, file_out2);
    file_out2.close();


    return 0;
}
// -----------------------------------------------------------
/*
Use this to plot (using gnuplot):

plot file using 1:2 title "Forward difference method", file using 1:3 title "Backward difference method", file using 1:4 title "Crank-Nicolson scheme"

plot file using 1:2 title "Forward difference method", file using 1:3 title "Backward difference method", file using 1:4 title "Crank-Nicolson scheme", file using 1:5 title "Exact Solution"

plot file using 1:2 title "Forward difference method" w lp, file using 1:3 title "Backward difference method" w lp, file using 1:4 title "Crank-Nicolson scheme" w lp, file using 1:5 title "Exact Solution"
*/
