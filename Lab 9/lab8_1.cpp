#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;
#define print(str, val)\
    cout << str << "\t:\t" << val << "\n";
#define epsilon -0.00001

// -----------------------------------------------------------
double f(double t, double y)
{
    return t*exp(3*t) - 2*y;
}

double orig_sol (double t)
{
    return t*exp(3*t)/5.0 - exp(3*t)/25.0 + exp(-2*t)/25.0;
}

double initial_value ()
{
    return 0.0;
}

double initial_t ()
{
    return 0.0;
}

double evaluate_t ()
{
    return 2.0;
}

double step_size ()
{
    return 0.2;
}
// -----------------------------------------------------------
double order_2_KH (double t_0, double t_1, double h)
{
    // y(t+1) = y(t) + h*f(y, t);
    double current_t = t_0;
    double y_0 = initial_value();
    double y_1;
    while (current_t < t_1 + epsilon)
    {
        y_1 = y_0 + h*f(current_t + 0.5*h, y_0 + 0.5 * h * f(current_t, y_0));
        y_0 = y_1;
        current_t += h;
    }
    return y_1;
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
int main()
{
    double t_0 = initial_t();
    double t_1 = evaluate_t();
    double h = step_size();

    vector<double> list_h;
    vector<double> list_error_order_2;
    vector<double> list_alpha_order_2;
    vector<double> list_error_order_4;
    vector<double> list_alpha_order_4;

    // For 7 steps of h from 0.1 to h/2 ... 7 times
    for (int i = 0; i < 7; i++)
    {
        // Order 2 KH Method at t = t_1
        double y_2 = order_2_KH (t_0, t_1, h);

        // Order 4 KH Method at t = t_1
        double y_4 = order_4_KH (t_0, t_1, h);

        // Calculate original solution
        double original_soln = orig_sol(evaluate_t());

        list_h.push_back (h);
        list_error_order_2.push_back (abs(y_2 - original_soln));
        list_error_order_4.push_back (abs(y_4 - original_soln));
        if (i > 0)
        {
            double alpha_2 = log(list_error_order_2[i-1]/list_error_order_2[i])/log(2.0);
            double alpha_4 = log(list_error_order_4[i-1]/list_error_order_4[i])/log(2.0);
            list_alpha_order_2.push_back (alpha_2);
            list_alpha_order_4.push_back (alpha_4);
        }
        h = h/2.0;
    }

    cout << "h\t\tError Order 2\tAlpha Order 2\tError Order 4\tAlpha Order 4\n";
    for (int i=0; i<7; i++)
    {
        cout << fixed << setprecision(10) << setw(10);
        cout << list_h[i] << "\t" << list_error_order_2[i] << "\t";
        if (i < list_alpha_order_2.size())
            cout << list_alpha_order_2[i] << "\t";
        else
            cout << "\t\t";
        cout << list_error_order_4[i] << "\t";
        if (i < list_alpha_order_4.size())
            cout << list_alpha_order_4[i] << "\n";
        else
            cout << "\t\t";
    }
    cout << "\n";
    return 0;
}
// -----------------------------------------------------------
