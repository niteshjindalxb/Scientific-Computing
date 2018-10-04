#include <iostream>
#include <cmath>
#include <vector>
#include <utility>

using namespace std;
#define print(str, val)\
    cout << str << "\t:\t" << val << "\n";

// -----------------------------------------------------------
// double f(double t, double y)
// {
//     return y/t - (y*y)/(t*t);
// }

// double orig_sol (double t)
// {
//     return (double)t/((double)1.0+ (double)log(t));
// }

// double initial_value ()
// {
//     return 1.0;
// }

// double initial_t ()
// {
//     return 1.0;
// }

// double evaluate_t ()
// {
//     return 2.0;
// }

// double step_size ()
// {
//     return 0.1;
// }
// -----------------------------------------------------------
// double f(double t, double y)
// {
//     return 1 + y/t + (y*y)/(t*t);
// }

// double orig_sol (double t)
// {
//     return t*tan(log(t));
// }

// double initial_value ()
// {
//     return 0.0;
// }

// double initial_t ()
// {
//     return 1.0;
// }

// double evaluate_t ()
// {
//     return 3.0;
// }

// double step_size ()
// {
//     return 0.2;
// }
// -----------------------------------------------------------
// double f(double t, double y)
// {
//     return -(y+1)*(y+3);
// }

// double orig_sol (double t)
// {
//     return -3.0 + 2.0/(1+exp(-2.0*t));
// }

// double initial_value ()
// {
//     return -2.0;
// }

// double initial_t ()
// {
//     return 0.0;
// }

// double evaluate_t ()
// {
//     return 2.0;
// }

// double step_size ()
// {
//     return 0.2;
// }
// -----------------------------------------------------------
double f(double t, double y)
{
    return -5*y + 5*t*t + 2*t;
}

double orig_sol (double t)
{
    return t*t + 1.0/3.0 * exp(-5.0*t);
}

double initial_value ()
{
    return 1.0/3.0;
}

double initial_t ()
{
    return 0.0;
}

double evaluate_t ()
{
    return 1.0;
}

double step_size ()
{
    return 0.1;
}
// -----------------------------------------------------------
double euler_approx (double t_0, double t_1, double h, vector<pair<double, double> > &my_pairs)
{
    // y(t+1) = y(t) + h*f(y, t);
    double current_t = t_0;
    double y_0 = initial_value();
    double y_1;
    while (current_t <= t_1)
    {
        my_pairs.push_back (make_pair(current_t, y_0));
        y_1 = y_0 + h*f(current_t, y_0);
        y_0 = y_1;
        current_t += h;
    }
    my_pairs.push_back (make_pair(current_t, y_0));
    return y_1;
}
// -----------------------------------------------------------
double linear_interpolation(pair<double, double> point_1, pair<double, double> point_2, double evaluate_at)
{
    return point_1.second + (evaluate_at - point_1.first) * (point_2.second - point_1.second) / (point_2.first - point_1.first);
}

// -----------------------------------------------------------
int main()
{
    double t_0 = initial_t();
    double t_1 = evaluate_t();
    double h = step_size();
    vector<pair<double, double> > my_pairs;

    // Euler approximation at t = t_0
    double y = euler_approx (t_0, t_1, h, my_pairs);
    print("Approximated Soln", y);
    double original_soln = orig_sol(evaluate_t());
    print("Original Soln\t", original_soln);

    // Now to find approximate solution at point_1, point_2 using linear interpolation
    double point[] = {1.05, 1.93};

    int count = 0;
    while (count < 2)
    {
        int pos = -1;
        // cout << "Points\n";
        for (int i = 0; i < my_pairs.size() - 1; i++)
        {
            // print(my_pairs[i].first, my_pairs[i].second);
            if (my_pairs[i+1].first > point[count])
            {
                pos = i;
                break;
            }
        }
        if (pos == -1)
        {
            cout << "Error Occurred\n";
            return 1;
        }
        cout << "------------\n";
        print (my_pairs[0].first, my_pairs[0].second);
        print (my_pairs[pos].first, my_pairs[pos].second);
        print (my_pairs[pos+1].first, my_pairs[pos+1].second);
        cout << "------------\n";
        double result = linear_interpolation (my_pairs[pos], my_pairs[pos+1], point[count]);
        print (point[count], result);
        print ("Real", orig_sol(point[count]));
        count ++;
    }
}
// -----------------------------------------------------------