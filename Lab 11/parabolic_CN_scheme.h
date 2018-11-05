#ifndef CN_SCHEME_H
#define CN_SCHEME_H true

std::vector<double> vector_stencil (std::vector<double> prev_soln, double lambda);
std::vector<std::vector<double> > CN_coef_matrix (double lambda);
std::vector<double> CN_method ();

#endif
