#ifndef MAIN_H
#define MAIN_H true
/* ********************************** */
#define print(str)\
    cout << str << "\n";
#define epsilon 1e-4
#define _USE_MATH_DEFINES
#define STABILITY 0.5

/* ********************************** */
const int x_grid_size = 5;
const int y_grid_size = 7;
const double x_range_from = 0.0;
const double x_range_to = 10.0;
const double y_range_from = 0.0;
const double y_range_to = 15.0;

double orig_sol (double x, double y);
double boundary_cond (double, double);
double initial_cond (double x);
std::vector<double> calc_exact_soln ();
std::vector<double> calc_error (std::vector<double> &approx_soln, std::vector<double> &exact_soln);

#endif
