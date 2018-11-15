#include <bits/stdc++.h>
#include "display.h"

void display (std::vector <double> &arr)
{
    cout << fixed << setprecision(10);
    for (int i = 0; i < arr.size(); i++)
        cout << arr[i] << "\n";
    cout << endl;
}
void display_2_vectors (std::vector <double> &arr1, std::vector <double> &arr2)
{
    cout << fixed << setprecision(10);
    for (int i = 0; i < arr1.size(); i++)
        cout << arr1[i] << "\t" << arr2[i] << endl;
}
void display_5_vectors (std::vector <double> &arr1, std::vector <double> &arr2, std::vector <double> &arr3, std::vector <double> &arr4, std::vector <double> &arr5, ofstream &file_out)
{
    file_out << fixed << setprecision(10);
    for (int i = 0; i < arr1.size(); i++)
        file_out << arr1[i] << "\t" << arr2[i] << "\t" << arr3[i] << "\t" << arr4[i] << "\t" << arr5[i] << endl;
}
void display_4_vectors (std::vector <double> &arr1, std::vector <double> &arr2, std::vector <double> &arr3, std::vector <double> &arr4, ofstream &file_out)
{
    file_out << fixed << setprecision(10);
    for (int i = 0; i < arr1.size(); i++)
        file_out << arr1[i] << "\t" << arr2[i] << "\t" << arr3[i] << "\t" << arr4[i] << endl;
}
