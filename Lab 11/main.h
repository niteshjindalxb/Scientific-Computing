#ifndef MAIN_H
#define MAIN_H true
/* ********************************** */
#define print(str)\
    cout << str << "\n";
#define epsilon 0.00000000001
#define _USE_MATH_DEFINES
#define STABILITY 0.5

/* ********************************** */
double orig_sol (double x, double t);
double x_step_size ();
double t_step_size ();
double x_range_from();
double x_range_to();
double t_range_from();
double evaluate_at();
double parabolic_c();
double boundary_cond ();
double initial_cond (double x);
std::vector<double> initial_value ();
std::vector<double> mesh_points ();
void display (std::vector <double> &arr);
void display_2_vectors (std::vector <double> &arr1, std::vector <double> &arr2);
void display_4_vectors (std::vector <double> &arr1, std::vector <double> &arr2, std::vector <double> &arr3, std::vector <double> &arr4, std::ofstream &file_out);
void display_5_vectors (std::vector <double> &arr1, std::vector <double> &arr2, std::vector <double> &arr3, std::vector <double> &arr4, std::vector <double> &arr5, std::ofstream &file_out);
std::vector<double> calc_exact_soln ();
std::vector<double> calc_error (std::vector<double> &approx_soln, std::vector<double> &exact_soln);

#endif
