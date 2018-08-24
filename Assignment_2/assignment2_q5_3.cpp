#include<iostream>
#include <iomanip>
#include <cmath>
using namespace std;
 
int main()
{
    double sum = 0;
    int n = 10000;

    for (int i = n; i >= 1; --i)
    {
    	double temp = 1.0/i;
    	temp = round(temp*10000);
    	temp = temp/10000.0;
        sum += temp;
    } 
    cout << setprecision (20);
    cout << "sum = " << sum << endl;
    return 0;
}
