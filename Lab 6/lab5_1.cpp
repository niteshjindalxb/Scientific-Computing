#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

#define tolerance 1e-5

double f (double x)
{
    return exp(2*x) * sin(3*x);
}

double simpson (double a, double b)
{
    return (b-a)*(f(a) + 4*f((a+b)/2) + f(b) ) / 6.0;
}

double adaptive_simpson (double a, double b, double max_error, vector<double> &x_points)
{
    // Calculate S1
    double integral_S1 = simpson (a, b);

    // Calculate S2
    double integral_S2 = simpson (a, (a+b)/2) + simpson ((a+b)/2, b);

    // Observed error
    double observed_error = abs(integral_S2 - integral_S1)/10.0;

    // If Observed error is less than max_error then return integral_S2
    if (observed_error <= max_error)
    {
        // Add x_points to vector
        x_points.push_back (a);
        x_points.push_back ((a+b)/2.0);
        x_points.push_back (b);
        return integral_S2;
    }

    vector<double> left_x_points;
    vector<double> right_x_points;

    // else, divide error in two sub-intervals and then proceed with adaptive simpson in each interval
    double left_adaptive_simpson = adaptive_simpson(a, (a+b)/2, max_error/2.0, left_x_points);
    double right_adaptive_simpson = adaptive_simpson((a+b)/2, b, max_error/2.0, right_x_points);

    for (int i=0; i<left_x_points.size() - 1; i++)
        x_points.push_back (left_x_points[i]);

    for (int i=0; i<right_x_points.size(); i++)
        x_points.push_back (right_x_points[i]);

    return left_adaptive_simpson + right_adaptive_simpson;
}

int main()
{
    double range_from = 1.0;
    double range_to = 3.0;

    vector<double> x_points;

    double integral = adaptive_simpson (range_from, range_to, tolerance, x_points);
    // cout << "Integral is = " << integral << endl;

    // cout << "x_size = " << x_points.size() << endl;
    for (int i=0; i<x_points.size(); i++)
    {
        cout << x_points[i] << "\t" << 0 << endl;
    }

    return 0;
}
