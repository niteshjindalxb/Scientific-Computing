#include <iostream>
#include <vector>
using namespace std;

// f(x)
double f(double x)
{
    return (double) 1.0/(1.0 + x*x);
}

// f'(x)
double df(double x)
{
    return (double) (-2 * x) / ((1.0 + x*x) * (1.0 + x*x));
}

// Lagrange Interpolation
vector<double> lagrange_poly (vector<double> &x, vector<double> &x_points)
{
    vector<double> y;
    for (int i=0; i<x.size(); i++)
    {
        double cur_y = 0.0;
        double numerator;
        double denominator;
        double f_x;

        // At how many points are you calculating your lagrange
        for (int j=0; j < x_points.size(); j++)
        {
            double cur_x = x[i];

            // point at which you're evaluating your lagrange
            double evaluate_x = x_points[j];

            numerator = 1.0;
            denominator = 1.0;
            f_x = f(evaluate_x);

            // All the terms in multiplication of numerator and denominator
            for (int k = 0; k < x_points.size(); k++)
            {
                double temp_x = x_points[k];
                if (temp_x == evaluate_x)
                    continue;

                numerator = numerator * (cur_x - temp_x);
                denominator = denominator * (evaluate_x - temp_x);
            }
            cur_y += numerator*f_x / denominator;
        }

        // push value of y corresponding to cur_y
        y.push_back (cur_y);
    }
    return y;
}

// vector<double> hermite_poly (vector<double> &x)
// {

// }

void display (vector<double> &x, vector<double> &y)
{
    for (int i=0; i<x.size(); i++)
        cout << x[i] << "\t" << y[i] << endl;
}

int main()
{
    vector<double> x_points;
    for (int i = -5; i <= 5; i++)
        x_points.push_back (i);

    vector<double> x;
    vector<double> y;
    for (int i = -5; i < 5; i++)
    {
        double h = 0.01;
        for (double j = 0; j < 100; j++)
            x.push_back (i + h*j);
    }
    x.push_back (5);

    y = lagrange_poly (x, x_points);

    display(x, y);
    return 0;
}
