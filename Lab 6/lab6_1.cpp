
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

#define MAX_ITERATION 53
#define print(x, y) \
	cout << "Iteration number = " << x << " and result = " << y << endl;

double f (vector<double> coef, double x)
{
    double result = coef[0];
	for (int i=0; i<coef.size(); i++)
        result = result*x + coef[i+1];
    return result;
}

double df (double x)
{
	return exp(x) - 1;
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

// void newton_method (double p0, int iteration_num, vector<double> &alpha)
// {
//     double denominator = df (p0);
//     double p1;

//     vector<double> calc_root;

//     while(iteration_num < MAX_ITERATION)
//     {
//         denominator = df (p0);
//         if (denominator == 0)
//         {
//             cout << "Slope is zero. Method stopped here\n";
//             calc_root.push_back(p0);
//             break;
//         }
//         p1 = p0 - f(p0)/(denominator);
//         // print(iteration_num, p1);
//         calc_root.push_back(p1);
//         if (abs(p1 - p0) < tolerance)
//             break;
//         else
//         {
//             p0 = p1;
//             iteration_num++;
//         }
//     }
//     for(int i=2; i<calc_root.size(); i++)
//     {
//         double alpha_value = log(calc_root[i-1]/calc_root[i])/log(calc_root[i-2]/calc_root[i-1]);
//         alpha.push_back (alpha_value);
//     }
// }

int main()
{
    double p0 = -1;
    int iteration_num = 1;
	// newton_method (p0, iteration_num, alpha);

    vector<double> coef = legendre_poly (5);
    for (int i=0; i<coef.size(); i++)
    {
        cout << coef[i] << " ";
    }
    cout << endl;

	return 0;
}