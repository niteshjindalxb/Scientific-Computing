#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <utility>
#include <iomanip>
using namespace std;

#define print(str, val)\
    cout << str << "\t:\t" << val << "\n";
#define epsilon -0.00000000000000000000001

int case_num = 4;
// -----------------------------------------------------------
double f(double t, double y)
{
    switch (case_num) {
        case 1: return y/t - (y*y)/(t*t);
        case 2: return 1.0 + y/t + (y*y)/(t*t);
        case 3: return -(y+1)*(y+3);
        case 4: return -5*y + 5*t*t + 2*t;
    }
}
double orig_sol (double t)
{
    switch (case_num) {
        case 1: return t/((double)1.0+ log(t));
        case 2: return t*tan(log(t));
        case 3: return -3.0 + 2.0/(1+exp(-2.0*t));
        case 4: return t*t + 1.0/3.0 * exp(-5.0*t);
    }
}
double initial_value ()
{
    switch (case_num) {
        case 1: return 1.0;
        case 2: return 0.0;
        case 3: return -2.0;
        case 4: return 1.0/3.0;
    }
}
double initial_t ()
{
    switch (case_num) {
        case 1: return 1.0;
        case 2: return 1.0;
        case 3: return 0.0;
        case 4: return 0.0;
    }
}
double evaluate_t ()
{
    switch (case_num) {
        case 1: return 2.0;
        case 2: return 3.0;
        case 3: return 2.0;
        case 4: return 1.0;
    }
}
double step_size ()
{
    switch (case_num) {
        case 1: return 0.1;
        case 2: return 0.2;
        case 3: return 0.2;
        case 4: return 0.1;
    }
}
// -----------------------------------------------------------
double euler_approx (double t_0, double t_1, double h, vector<pair<double, double> > &my_pairs, bool store_pairs)
{
    // y(t+1) = y(t) + h*f(y, t);
    double current_t = t_0;
    double y_0 = initial_value();
    double y_1;
    while (current_t < t_1 + epsilon)
    {
        if (store_pairs)
            my_pairs.push_back (make_pair(current_t, y_0));
        y_1 = y_0 + h*f(current_t, y_0);
        y_0 = y_1;
        current_t += h;
    }
    if (store_pairs)
        my_pairs.push_back (make_pair(current_t, y_0));
    return y_1;
}
// -----------------------------------------------------------
double linear_interpolation(pair<double, double> point_1, pair<double, double> point_2, double evaluate_at)
{
    return point_1.second + (evaluate_at - point_1.first) * (point_2.second - point_1.second) / (point_2.first - point_1.first);

}
// -----------------------------------------------------------
void display(std::vector<double> &list_h, std::vector<double> &list_error, std::vector<double> &list_alpha)
{
    cout << "h\t\tError\t\tAlpha\n";
    for (int i = 0; i < list_h.size(); i++)
    {
        cout << fixed << setprecision(10) << setw(10);
        cout << list_h[i] << "\t" << list_error[i] << "\t";
        if (i < list_alpha.size())
            cout << list_alpha[i] << "\n";
        else
            cout << "\n";
    }
    cout << "\n";
}
// -----------------------------------------------------------
void calc_order_of_convergence (std::vector<double> &list_h, std::vector<double> &list_error, std::vector<double> &list_alpha)
{
    list_h.clear();
    list_error.clear();
    list_alpha.clear();
    vector<pair<double, double> > _;
    int max_iterations = 20;
    double h = step_size();
    for (int i = 0; i < max_iterations; i++)
    {
        // Euler approximation method at t = t_1
        double y = euler_approx(initial_t(), evaluate_t(), h, _, /*store_pairs = */ false);

        // Calculate original solution
        double original_soln = orig_sol(evaluate_t());

        list_h.push_back (h);
        list_error.push_back (abs(y - original_soln));
        if (i > 0)
        {
            double alpha = log(list_error[i-1]/list_error[i])/log(2.0);
            list_alpha.push_back (alpha);
        }
        h = h/2.0;
    }
}
// -----------------------------------------------------------
void write_to (string filename, std::vector<double> &arr1, std::vector<double> &arr2)
{
    std::ofstream file(filename, std::ios::out);
    if (file.is_open()) {
        for (size_t i = 0; i < max(arr1.size(), arr2.size()); i++) {
            if (i < arr1.size())
                file << arr1[i] << "\t";
            else
                file << "\t";

            if (i < arr2.size())
                file << arr2[i] << "\n";
            else
                file << "\n";
        }
        file.close();
    }
}
// -----------------------------------------------------------
int main()
{
    double t_0 = initial_t();
    double t_1 = evaluate_t();
    double h = step_size();
    vector<pair<double, double> > my_pairs;

    /*********** Find Euler approximation at t = t_0 ***********/
    double y = euler_approx (t_0, t_1, h, my_pairs, /*store_pairs = */ true);
    print("Approximated Soln", y);
    double original_soln = orig_sol(evaluate_t());
    print("Original Soln\t", original_soln);

    /*********** Find approximate order of convergence ***********/
    // For 7 steps of h from h0 to h0/2 ... 7 times
    vector<double> list_h;
    vector<double> list_error;
    vector<double> list_alpha;

    calc_order_of_convergence(list_h, list_error, list_alpha);
    string filename;
    switch (case_num) {
        case 1: filename = "caseA.txt";
                break;
        case 2: filename = "caseB.txt";
                break;
        case 3: filename = "caseC.txt";
                break;
        case 4: filename = "caseD.txt";
                break;
    }
    write_to (filename, list_h, list_error);
    display (list_h, list_error, list_alpha);

    /*********** Find approximate solution at point_1, point_2 using linear interpolation ***********/
    double point[] = {
                        1.05, 1.93,
                        2.1 , 2.75,
                        1.3 , 1.93,
                        0.54, 0.94
                    };

    int count = 0;
    while (count < 2)
    {
        int pos = -1;
        for (int i = 0; i < my_pairs.size() - 1; i++)
        {
            if (my_pairs[i+1].first > point[count + (case_num - 1)*2] - epsilon)
            {
                pos = i;
                break;
            }
        }
        if (pos == -1)
        {
            cout << "Error Occurred\n";
            count ++;
            continue;
        }
        double result = linear_interpolation (my_pairs[pos], my_pairs[pos+1], point[count + (case_num - 1)*2]);
        print (point[count + (case_num - 1)*2], result);
        print ("Real\t", orig_sol(point[count + (case_num - 1)*2]));
        count ++;
    }
    return 0;
}
// -----------------------------------------------------------
