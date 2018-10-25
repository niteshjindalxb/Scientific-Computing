
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
// Display Polynomial
void display (vector<double> &arr)
{
    cout << "Polynomial :\n";
    for (size_t i = 0; i < arr.size(); i++)
        cout << arr[i] << " ";
    cout << "\n";
}

int main()
{
    double p0 = -1;
    int iteration_num = 1;
    int degree = 2;

    vector<double> coef = legendre_poly (degree);
    display(coef);

	return 0;
}