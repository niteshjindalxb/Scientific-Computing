/*
Question 1-
Let f(x) = e^x − x − 1 . Use Newton’s method to find the zero in [−1,1]. Compare the results with those obtained using Secant method and bisection method. In all the cases  compute the root up to an accuracy of 10^−6. */

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
using namespace std;

#define tolerance 0.000001
#define print(x, y) \
	cout << "Iteration number = " << x << " and result = " << y << endl;

double f (double x)
{
	return exp(x)-x-1;
}

double df (double x)
{
	return exp(x)-1;
}

double newton_method (double old_result, int iteration_num, vector<double> &file)
{
	double denominator = df (old_result);
	if (denominator == 0)
	{
		cout << "Slope is zero. Method stopped here\n";
		return old_result;
	}
	double new_result = old_result - f(old_result)/(denominator);
	// print(iteration_num, new_result);
	file.push_back(new_result);
	if (abs(new_result - old_result) < tolerance)
		return new_result;
	else
		return newton_method (new_result, iteration_num+1, file);
}

double bisection_method (double range_from, double range_to, int iteration_num, vector<double> &file)
{
	if (range_to - range_from == 0)
		return range_from;
	
	double mid = range_from + (range_to-range_from)/2.0;
	double temp = df(mid);

	// print(iteration_num, mid);
	file.push_back(mid);
	if (temp == 0)
		return mid;
	else if (temp < 0)
		return bisection_method (mid, range_to, iteration_num+1, file);
	else
		return bisection_method (range_from, mid, iteration_num+1, file);
}

double secant_method (double x_n_1, double x_n_2, int iteration_num, vector<double> &file)
{
	double x_n;
	while (abs(x_n_1 - x_n_2) > tolerance)
	{
		if (f(x_n_1) == f(x_n_2))
		{
			cout << "Slope is zero. Method stopped here\n";
			break;
		}
		x_n = x_n_1 - f(x_n_1)* (x_n_2-x_n_1)/(f(x_n_2)-f(x_n_1));
		file.push_back(x_n);
		// print (iteration_num++, x_n);
		x_n_2 = x_n_1;
		x_n_1 = x_n;
	}
	return x_n;
}

int main()
{
	vector<double> newton_result;
	vector<double> bisection_result;
	vector<double> secant_result;

	double range_from = -1.0;
	double range_to = 1.0;

	// cout << "\nSecant Method starts\n";
	double root = secant_method (range_from, range_to, 1, secant_result);
	// cout << "Root = " << root << endl;

	// cout << "Newton Method starts\n";
	root = newton_method (0.1, 1, newton_result);
	// cout << "Root = " << root << endl;

	// cout << "\nBisection Method starts\n";
	root = bisection_method (range_from, range_to, 1, bisection_result);
	// cout << "Root = " << root << endl;
	
	cout << "Iteration number \t Newton Method \t\t Secant Method \t\t Bisection Method\n";
	int size1 = newton_result.size();
	int size2 = secant_result.size();
	int size3 = bisection_result.size();
	int max_size = max (size1, size2);
	max_size = max (max_size, size3);
	
	// cout << size1 << endl;
	// cout << size2 << endl;
	// cout << size3 << endl << endl;
    cout << fixed << setprecision(6);
	for (int i=0; i<max_size; i++)
	{
		cout << i+1 << "\t\t\t";
		if (size1 > i)
			cout << newton_result[i] << "\t\t";
		else
			cout << "\t" << "\t\t";

		if (size2 > i)
			cout << secant_result[i] << "\t\t";
		else
			cout << "\t" << "\t\t";

		if (size3 > i)
			cout << bisection_result[i] << "\t\t";
		else
			cout << "\t" << "\t\t";
		
		cout << endl;
	}

	return 0;
}
