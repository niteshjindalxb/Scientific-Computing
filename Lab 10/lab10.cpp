#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;
#define print(str)\
    cout << str << "\n";
#define epsilon -0.00000000001
#define tolerance 0.000001
#define MAX_ITERATION 53
#define _USE_MATH_DEFINES

// -----------------------------------------------------------
double f (double t, double y)
{
    return exp(y);
}
double step_size ()
{
    return 0.01;
}
double g(double w_i, double w_i_1, double w_i_2, double w)
{
    return w_i + step_size()/24.0 * (9.0*exp(w) + 19*exp(w_i) - 5*exp(w_i_1) + exp(w_i_2));
}
double orig_sol (double t)
{
    return 1.0 - log(1 - M_E*t);
}
double error (double pn, double p)
{
    return abs(pn-p);
}
double initial_t ()
{
    return 0.0;
}
double evaluate_at()
{
    return 0.2;
}
double initial_value ()
{
    return 1.0;
}
// -----------------------------------------------------------
void display (vector <double> &arr)
{
    cout << fixed << setprecision(10);
    for (int i = 0; i < arr.size(); i++)
        cout << arr[i] << "\n";
    cout << endl;
}
// -----------------------------------------------------------
void display_2_vectors (vector <double> &arr1, vector <double> &arr2)
{
    cout << fixed << setprecision(10);
    for (int i = 0; i < arr1.size(); i++)
        cout << arr1[i] << "\t" << arr2[i] << endl;
}
// -----------------------------------------------------------
// Computes next w value and append it at last of w_values list
void fixed_point_method (vector <double> &w_values, vector <double> &fixed_alpha)
{
    // To store inner alpha values
    vector <double> iter_alpha;

    // To calculate w[i+1] using fixed point iteration
    double p1;
    // Take starting point as w[i]
    double p0 = w_values [w_values.size()-1];

    // To calculate alpha values
    vector <double> fixed_point_value;

    // g requires 3 previous points to evaluate w
    double w_i = w_values[w_values.size()-1];
    double w_i_1 = w_values[w_values.size()-2];
    double w_i_2 = w_values[w_values.size()-3];

    // Fixed point method starts
    int iteration_num = 0;
    while(iteration_num < MAX_ITERATION)
    {
        p1 = g(w_i, w_i_1, w_i_2, p0);
        fixed_point_value.push_back(p1);
        if (abs(p1 - p0) < tolerance)
            break;
        else
        {
            p0 = p1;
            iteration_num++;
        }
    }
    // We got the fixed point i.e. p1
    w_values.push_back (p1);
    
    // Now calculate alpha values
    double original_soln = orig_sol (step_size() * w_values.size());
    for(int i = 3; i < fixed_point_value.size(); i++)
    {
        double alpha_value = log(error(fixed_point_value[i-1], original_soln)/error(fixed_point_value[i], original_soln))/log(error(fixed_point_value[i-2], original_soln)/error(fixed_point_value[i-1], original_soln));
        iter_alpha.push_back (alpha_value);
    }
    double sum_alpha = 0.0;
    for (int i=0; i<iter_alpha.size(); i++)
        sum_alpha += iter_alpha[i];
    
    fixed_alpha.push_back (sum_alpha/iter_alpha.size());
}
// -----------------------------------------------------------
double order_4_KH (double t_0, double t_1, double h)
{
    // y(t+1) = y(t) + h*f(y, t);
    double current_t = t_0;
    double y_0 = initial_value();
    double y_1;
    while (current_t < t_1 + epsilon)
    {
        double k1 = f(current_t, y_0);
        double k2 = f(current_t+0.5*h, y_0 + 0.5* h* k1);
        double k3 = f(current_t+0.5*h, y_0 + 0.5* h* k2);
        double k4 = f(current_t+h, y_0 + h* k3);
        y_1 = y_0 + (double)1.0/6.0 * h * (k1 + 2*k2 + 2*k3 + k4);
        y_0 = y_1;
        current_t += h;
    }
    return y_1;
}
// -----------------------------------------------------------
// Computes next w value and append it at last of w_values list
void adams_PC (vector <double> &w_values)
{
    double h = step_size();
    // Take last 4 values and compute next 5th value
    int size = w_values.size();
    double val0 = w_values[size-4];
    double val1 = w_values[size-3];
    double val2 = w_values[size-2];
    double val3 = w_values[size-1];

    double t0 = (size-4)*h;
    double t1 = (size-3)*h;
    double t2 = (size-2)*h;
    double t3 = (size-1)*h;
    double t4 = size*h;

    double val4_temp = val3 + h/24.0 * (55.0*f(t3, val3) - 59.0*f(t2, val2) + 37.0*f(t1, val1) - 9.0*f(t0, val0));
    double val4 = val3 + h/24.0 * (9*f(t4, val4_temp) + 19.0*f(t3, val3) - 5.0*f(t2, val2) + f(t1, val1));

    w_values.push_back (val4);
}
// -----------------------------------------------------------
void adams_moulton_implicit (vector <double> &w_values, vector <double> &fixed_alpha)
{
    double h = step_size();
    // For 20 steps of h from 0.0 to 0.0+h ... 20 times
    for (int i = 0; i < 20; i++)
    {
        if (i < 3)
            w_values.push_back (orig_sol(h*i));
        else
            fixed_point_method(w_values, fixed_alpha);
    }
}
// -----------------------------------------------------------
void calc_exact_soln (vector <double> &exact_soln)
{
    double h = step_size();
    for (int i = 0; i < 20; i++)
        exact_soln.push_back (orig_sol(h*i));
}
// -----------------------------------------------------------
void adams_predictor_corrector (vector <double> &w_values)
{
    double h = step_size();
    // Now to calculate approximate value of differential equation y' = e^y at t = 0.2
    // using Adams fourth-order predictor-corrector method
    // Firstly calculating 4 initial points using 4th order KH Method
    w_values.push_back(1.0);
    for (int i = 1; i < 4; i++)
        w_values.push_back(order_4_KH(initial_t(), h*i, h));

    // Now calculate remaining values using Adams fourth order predictor-corrector method
    double evaluate_at = 0.2;
    double current_at = 3*h;

    // Predictor-corrector method
    while (current_at < evaluate_at + epsilon)
    {
        adams_PC (w_values);
        current_at = h * (w_values.size()-1);
    }
}
// -----------------------------------------------------------
vector<double> calc_error (vector<double> &approx_soln, vector<double> &exact_soln)
{
    vector<double> list_error;
    for (int i = 0; i < exact_soln.size(); i++)
        list_error.push_back (error(approx_soln[i], exact_soln[i]));
    return list_error;
}
// -----------------------------------------------------------
int main()
{
    double t_0 = initial_t();
    double h = step_size();

    // Adams moulton implicit method
    vector <double> w_values_implicit;
    vector <double> fixed_alpha;
    // Exact Solution
    vector <double> exact_soln;
    vector <double> w_values_PC;

    adams_moulton_implicit (w_values_implicit, fixed_alpha);
    // Now to calculate approximate value of differential equation y' = e^y at t = 0.2
    // using Adams fourth-order predictor-corrector method
    adams_predictor_corrector (w_values_PC);
    // Calculate exact solution
    calc_exact_soln (exact_soln);

    // Error calculation
    vector <double> adams_moulton_implicit_error = calc_error (w_values_implicit, exact_soln);
    vector <double> w_values_PC_error = calc_error (w_values_PC, exact_soln);
    
    print("---------------------------------");
    print("Fixed Point\tPC");
    print("---------------------------------");
    display_2_vectors (adams_moulton_implicit_error, w_values_PC_error);
    print("---------------------------------");

    // Now to calculate order of convergence of PC Method
    vector <double> list_alpha_implicit;
    vector <double> list_alpha_PC;
    vector <double> list_error_implicit;
    vector <double> list_error_PC;
    int max_iteration = 10;
    for (int i = 0; i < max_iteration; i++)
    {
        w_values_implicit.clear();
        fixed_alpha.clear();
        w_values_PC.clear();

        adams_moulton_implicit (w_values_implicit, fixed_alpha);
        adams_predictor_corrector (w_values_PC);
        int size = w_values_implicit.size();
        int original_soln = orig_sol(evaluate_at());

        list_error_implicit.push_back (error(w_values_implicit[size-1], original_soln));
        list_error_PC.push_back (error(w_values_PC[size-1], original_soln));

        if (i > 0)
        {
            double temp = list_error_implicit[i-1]/list_error_implicit[i];
            list_alpha_implicit.push_back (log2(temp));

            temp = list_error_PC[i-1]/list_error_PC[i];
            list_alpha_PC.push_back (log2(temp));
        }
        h = h/2.0;
    }

    display_2_vectors (list_alpha_implicit, list_alpha_PC);    
    return 0;
}
// -----------------------------------------------------------