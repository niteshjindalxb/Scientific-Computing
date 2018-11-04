#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

// f(x)
double f(double x)
{
    return (double) 1.0/(1.0 + x*x);
}


// Lagrange Interpolation
vector<double> lagrange_poly (vector<double> &x, vector<double> &x_points, vector<double> &y_points)
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
            f_x = y_points[j];

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


void display (vector<double> &x, vector<double> &y)
{
    for (int i=0; i<x.size(); i++)
        cout << fixed << setprecision(10) << x[i] << "\t" << y[i] << endl;
}

int main()
{
    vector<double> x_points;
    x_points.push_back (0.3);
    x_points.push_back (0.4);
    x_points.push_back (0.5);
    x_points.push_back (0.6);

    vector<double> e_x;
    e_x.push_back (0.740818);
    e_x.push_back (0.670320);
    e_x.push_back (0.606531);
    e_x.push_back (0.548812);

    vector<double> y_points(x_points.size());
    for (int i=0; i<4; i++)
        y_points[i] = x_points[i] - e_x[i];

    vector<double> x;
    vector<double> y;
    x.push_back (0);

    y = lagrange_poly (x, y_points, x_points);
    display(x, y);
    return 0;
}
