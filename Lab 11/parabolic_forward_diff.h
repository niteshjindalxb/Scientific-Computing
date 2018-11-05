#ifndef FORWARD
#define FORWARD true

double forward_stencil (double x, double t, std::vector<double> prev_soln);
std::vector<double> forward_diff_method ();

#endif
