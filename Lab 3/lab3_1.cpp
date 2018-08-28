/*
Question 1-
Let f(x) = e^x − x − 1. Compute the approximate order of convergence of both Newton’s and Modified Newton’s
methods in finding a zero of f in the interval [− 1 , 1]. */

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
using namespace std;

#define tolerance 0.0001
#define MAX_ITERATION 53
#define _USE_MATH_DEFINES
#define print(x, y) \
	cout << "Iteration number = " << x << " and result = " << y << endl;
#define ROOT 0

double f (double x)
{
	return exp(x) - x - 1;
}

double df (double x)
{
	return exp(x) - 1;
}

double ddf (double x)
{
	return exp(x);
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
			std::cout << calc_root[i] << '\n';
		}
		std::cout << "Alpha_Value" << '\n';
    for(int i=2; i<calc_root.size(); i++)
    {
        double alpha_value = log(calc_root[i-1]/calc_root[i])/log(calc_root[i-2]/calc_root[i-1]);
        alpha.push_back (alpha_value);
    }
}

double modified_newton_function (double x)
{
    return x - (f(x) * df(x)) / ( df(x)*df(x) - f(x)*ddf(x) );
}

void modified_newton_method (double p0, int iteration_num, vector<double> &alpha)
{
    alpha.clear();
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
        p1 = modified_newton_function (p0);
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
			std::cout << calc_root[i] << '\n';
		}
		std::cout << "Alpha_Value" << '\n';
    for(int i=2; i<calc_root.size(); i++)
    {
        double alpha_value = log(calc_root[i-1]/calc_root[i])/log(calc_root[i-2]/calc_root[i-1]);
        alpha.push_back (alpha_value);
    }
}


int main()
{
    vector<double> alpha, modified_alpha;

	cout << "Newton Method starts\n";
    double p0 = -1;
    int iteration_num = 1;
	newton_method (p0, iteration_num, alpha);
    // cout << "alpha for newton's method\n";
    for(int i=0; i<alpha.size(); i++)
        cout << alpha[i] << endl;

	std::cout << "modified_newton_method" << '\n';

	modified_newton_method (p0, iteration_num, modified_alpha);
    // cout << "modified_alpha for modified newton's method\n";
    for(int i=0; i<modified_alpha.size(); i++)
        cout << modified_alpha[i] << endl;

	return 0;
}
