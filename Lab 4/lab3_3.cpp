/*
Question 3-
Let f(x) = (x-1)*(x-6)*(x-8), compare the approximate order of convergence of both secant and Newton’s method
in finding the root of f in the interval [0, 2] */

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
using namespace std;

#define tolerance 0.0001
#define MAX_ITERATION 53
#define _USE_MATH_DEFINES
#define ROOT 1

double f (double x)
{
	return (x-1)*(x-6)*(x-8);
}

double df (double x)
{
	return (x-1)*(x-6) + (x-6)*(x-8) + (x-1)*(x-8);
}

double error (double pn)
{
    return abs(pn-ROOT);
}

void newton_method (double p0, int iteration_num, vector<double> &alpha)
{
    double denominator = df (p0);
    double p1;

    vector<double> calc_root;

    while(iteration_num < MAX_ITERATION)
    {
        denominator = df (p0);
        if (denominator == 0)
        {
            cout << "Slope is zero. Method stopped here\n";
            calc_root.push_back(p0);
            break;
        }
        p1 = p0 - f(p0)/(denominator);
        // print(iteration_num, p1);
        calc_root.push_back(p1);
        if (abs(p1 - p0) < tolerance)
            break;
        else
        {
            p0 = p1;
            iteration_num++;
        }
    }
		for (size_t i = 0; i < calc_root.size(); i++) {
			std::cout << error(calc_root[i]) << '\n';
		}
		std::cout << "Alpha Values" << '\n';
    for(int i=2; i<calc_root.size(); i++)
    {
        double alpha_value = log(error(calc_root[i-1])/error(calc_root[i]))/log(error(calc_root[i-2])/error(calc_root[i-1]));
        alpha.push_back (alpha_value);
    }
}

double secant_method (double x_n_1, double x_n_2, int iteration_num, vector<double> &alpha)
{
	double x_n;
    vector<double> calc_root;

	while (abs(x_n_1 - x_n_2) > tolerance)
	{
		if (f(x_n_1) == f(x_n_2))
		{
			cout << "Slope is zero. Method stopped here\n";
			break;
		}
		x_n = x_n_1 - f(x_n_1)* (x_n_2-x_n_1)/(f(x_n_2)-f(x_n_1));
		calc_root.push_back(x_n);
		// print (iteration_num++, x_n);
		x_n_2 = x_n_1;
		x_n_1 = x_n;
	}
	for (size_t i = 0; i < calc_root.size(); i++) {
		std::cout << error(calc_root[i]) << '\n';
	}
	std::cout << "Alpha Values" << '\n';
    for(int i=2; i<calc_root.size(); i++)
    {
        double alpha_value = log(error(calc_root[i-1])/error(calc_root[i]))/log(error(calc_root[i-2])/error(calc_root[i-1]));
        alpha.push_back (alpha_value);
    }
	return x_n;
}

int main()
{
    vector<double> newton_alpha, secant_alpha;

	cout << "Newton Method starts\n";
    double p0 = -1;
    int iteration_num = 1;
	newton_method (p0, iteration_num, newton_alpha);
    // cout << "newton_alpha for newton's method\n";
    for(int i=0; i<newton_alpha.size(); i++)
        cout << newton_alpha[i] << endl;

    double x_n_1 = 0.5;
    double x_n_2 = 1.5;

	std::cout << "secant_method" << '\n';
	secant_method (x_n_1, x_n_2, iteration_num, secant_alpha);
    cout << "secant_alpha for secant's method\n";
    for(int i=0; i<secant_alpha.size(); i++)
        cout << secant_alpha[i] << endl;

	return 0;
}
