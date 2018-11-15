#include <bits/stdc++.h>
#include "gauss_seidel.h"
using namespace std;

vector<double> initial_soln(int order_matrix)
{
    std::vector<double> init_soln(order_matrix, 0.0);
    return init_soln;
}

std::vector<double> gauss_seidel(vector<vector<double> > A, vector<double> b)
{
    int size = b.size();
	// Declare and Initialize starting solution vector
    std::vector<double> x = initial_soln(b.size());
	double error = -1.0;

	// Keep count of iteration
	int iter = 0;
	do
	{
		iter++;
		for (int i = 0; i < A.size(); ++i)
		{
			double sum = 0.0;
			for (int j = 0; j < A.size(); ++j)
				if(i != j)
					sum -= A[i][j] * x[j];
			sum += b[i];
			double temp_soln = sum/A[i][i];

			// Calculate new error
			error = max(error, abs(x[i]-temp_soln));

			x[i] = temp_soln;
		}
		if(iter > MAX_ITERATIONS)
			break;
	} while(error > tolerance);
	return x;
}
