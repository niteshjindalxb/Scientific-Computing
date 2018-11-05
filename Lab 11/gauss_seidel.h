#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H true

#define MAX_ITERATIONS 1000
#define tolerance 0.001

std::vector<double> initial_soln(int order_matrix);
std::vector<double> gauss_seidel(std::vector<std::vector<double> > A, std::vector<double> b);

#endif
