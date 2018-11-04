#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>
using namespace std;

void display (vector<double> &arr)
{
    for (int i = 0; i < arr.size(); i++)
        cout << arr[i] << " ";
    cout << "\n";
}

double orig_soln (double x)
{
    return 0.5118277;
}

double evaluate_at ()
{
    return 1.5;
}
double forward_product_terms (double x, vector<pair<double, double> > &points, int upper)
{
    double result = 1.0;
    for (int i = 0; i <= upper; i++)
        result = result * (x - points[i].first);
    return result;
}
double backward_product_terms (double x, vector<pair<double, double> > &points, int lower)
{
    double result = 1.0;
    for (int i = lower; i <= points.size()-1; i++)
        result = result * (x - points[i].first);
    return result;
}
double newton_forward_diff_method (vector<pair<double, double> > &points, double evaluate_at)
{
    vector <double> arr(points.size());
    double result = 0.0;

    for (int j = 0; j < points.size(); j++)
    {
        for (int i = points.size() - 1; i >= j; i--)
        {
            // Computing f[xi, .., xj]
            if (j == 0)
                arr [i] = points[i].second;
            else
                arr [i] = (arr[i] - arr[i-1])/(points[i].first - points[i-j].first);
        }
        // f[x0, x1, .., xj] = arr[j]
        result += arr[j] * forward_product_terms (evaluate_at, points, /* x0, .., xj-1 */ j-1);
    }
    return result;
}
double newton_backward_diff_method (vector<pair<double, double> > &points, double evaluate_at)
{
    vector <double> arr(points.size());
    double result = 0.0;

    for (int j = 0; j < points.size(); j++)
    {
        for (int i = points.size()-1; i >= j; i--)
        {
            // Computing f[xi, .., xj]
            if (j == 0)
                arr [i] = points[i].second;
            else
                arr [i] = (arr[i] - arr[i-1])/(points[i].first - points[i-j].first);
        }
        // f[x0, x1, .., xj] = arr[j]
        int n = points.size();
        result += arr[n-1] * backward_product_terms (evaluate_at, points, /* x(n-j), .., xn */ n-j);
    }
    return result;
}
int main()
{
    // map <double, double> f;
    vector<pair<double, double> > points;

    points.push_back (make_pair(1.0, 0.7651977));
    points.push_back (make_pair(1.3, 0.6200860));
    points.push_back (make_pair(1.6, 0.4554022));
    points.push_back (make_pair(1.9, 0.2818186));
    points.push_back (make_pair(2.2, 0.1103623));
    points.push_back (make_pair(2.5, 0.0483838));

    double forward_result = newton_forward_diff_method (points, evaluate_at());
    cout << fixed << setprecision(10) << forward_result << endl;

    double backward_result = newton_backward_diff_method (points, evaluate_at());
    cout << backward_result << endl;

    return 0;
}
