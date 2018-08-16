#include<iostream>
using namespace std;
 
int main()
{
    double sum = 0;
    int n = 10000;

    for (int i = 1; i <= n; ++i)
    {
        sum += 1.0/i;
    } 
    cout << "sum = " << sum << endl;
    return 0;
}
