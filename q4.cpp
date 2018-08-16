#include <iostream>
using namespace std;

int main(int argc, char const *argv[])
{
	double e = 1.0/2.0;
	double e_new = 1.0/2.0;
	while ((1+e_new) > 1)
	{
		e = e_new;
		e_new = e_new / 2.0;
	}
	cout << "e = " << e << endl;
	return 0;
}
