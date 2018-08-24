/*
Question 2-
Compute the approximate order of convergence of the fixed point iteration in finding the fixed point of the function cos(x) in [0 ,π/2]. */ 

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
using namespace std;

#define tolerance 0.000001
#define MAX_ITERATION 50
#define _USE_MATH_DEFINES

double f (double x)
{
	return cos(x);
}

double error (double pn, double p)
{
    return abs(pn-p);
}

void fixed_point_method (double p0, int iteration_num, vector<double> &alpha)
{
    double p1;
    vector<double> calc_root;

    while(iteration_num < MAX_ITERATION)
    {
        p1 = f(p0);
        calc_root.push_back(p1);
        if (abs(p1 - p0) < tolerance)
            break;
        else
        {
            p0 = p1;
            iteration_num++;
        }
    }
    for(int i=3; i<calc_root.size(); i++)
    {
        double alpha_value = log(error(calc_root[i-1], calc_root[i-2])/error(calc_root[i], calc_root[i-1]))/log(error(calc_root[i-2], calc_root[i-3])/error(calc_root[i-1], calc_root[i-2]));
        alpha.push_back (alpha_value);
    }
}

int main()
{
    vector<double> alpha;

	cout << "Fixed point Method starts\n";
    double p0 = M_PI/4;
    int iteration_num = 1;
	fixed_point_method (p0, iteration_num, alpha);
    for(int i=0; i<alpha.size(); i++)
        cout << alpha[i] << endl;
	return 0;
}
