#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
using namespace std;

#define tolerance 0.000000001
#define epsilon 0.0000001
#define print(val) \
	cout << fixed << setprecision(10) << val << endl;


double f (vector<double> coef, double x)
{
    double result = coef[0];
	for (int i=0; i<coef.size() - 1; i++)
        result = result*x + coef[i+1];
    return result;
}
double bisection_method (vector<double> coef, double range_from, double range_to)
{
	if (range_to - range_from < tolerance)
		return range_from;

	double mid = range_from + (range_to-range_from)/2.0;
	double mid_value = f(coef, mid);

	if (mid_value < tolerance && mid_value > -tolerance)
		return mid;
	else if (mid_value < 0)
		return bisection_method (coef, mid, range_to);
	else
		return bisection_method (coef, range_from, mid);
}
// Display Polynomial
void display (vector<double> &coef)
{
	print("--------------------------");
    for (size_t i = 0; i < coef.size(); i++)
        cout << coef[i] << " ";
    cout << "\n";
}
vector<double> legendre_poly (int deg)
{
    vector<double> coef;
    if (deg == 0)
    {
        coef.push_back (1.0);
        return coef;
    }
    else if (deg == 1)
    {
        coef.push_back (1.0);
        coef.push_back (0.0);
        return coef;
    }
    else
    {
        int j = deg-1;
        vector<double> coef_j = legendre_poly (j);
        vector<double> coef_j_1 = legendre_poly (j - 1);

        for (int i = 0; i < coef_j.size(); i++)
            coef_j[i] = ((double)(2*j + 1)/(double)(j + 1))*coef_j[i];

        coef_j.push_back (0.0);

        for (int i = 0; i < coef_j_1.size(); i++)
            coef_j_1[i] = ((double)j/(double)(j + 1))*coef_j_1[i];

        // Subtract two polynomials
        int k = 0;
        for (int i = 0; i < coef_j.size(); i++)
        {
            if (i < coef_j.size() - coef_j_1.size())
                continue;
            coef_j[i] = coef_j[i] - coef_j_1[k++];
        }
        return coef_j;
    }
}
set<double> find_roots (std::vector<double> coef)
{
	// We know that all of its roots occurs between (-1, 1)
	set<double> roots;
	int count = 0;
	double range_from = -1.0;
	double h = 2.0;

	while (count < coef.size() - 1 && h > 0.00001)
	{
		range_from = -1.0;
		for (double range_to = -1.0; range_to < 1.0 + epsilon; range_to += h)
		{
			double left = f (coef, range_from);
			double right = f (coef, range_to);
			if (left*right < 0)
			// Root exists in this interval
			{
				double new_root = bisection_method (coef, range_from, range_to);
				if (roots.find(new_root) == roots.end())
				{
					roots.insert (new_root);
					count ++;
				}
			}
			range_from = range_to;
		}
		h /= 2.0;
	}
	return roots;
}

int main()
{
	// Set degree of polynomial here
	int degree = 10;

	// Find coefficient of the polynomial
	vector<double> coef = legendre_poly (degree);

	// Find the roots of legendre polynomial
	set <double> set_root = find_roots (coef);
	vector<double> roots(set_root.begin(), set_root.end());

	// Display roots
	display(roots);

	return 0;
}
