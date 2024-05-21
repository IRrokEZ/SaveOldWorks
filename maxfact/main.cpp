#include <iostream>
#include <limits>
#include <iomanip>

using namespace std;

int main()
{
    //int64_t fact, i;
    //int64_t res;
    //double a;
    //long double b;
    //long long double c;
    //cout << endl << sizeof(a) << " " << sizeof(b) << endl; //" " << sizeof(c) << endl;
    //return 0;
    //cout << endl << abs(-5) << endl;
    //return 0;
    long double t = 123.123456789;
    cout << setprecision(0) << t << endl;
    cout << setprecision(1) << t << endl;
    cout << setprecision(2) << t << endl;
    cout << setprecision(3) << t << endl;
    cout << setprecision(4) << t << endl;
    cout << setprecision(5) << t << endl;
    cout << setprecision(6) << t << endl;
    cout << setprecision(7) << t << endl;
    cout << setprecision(8) << t << endl;
    cout << setprecision(9) << t << endl;
    return 0;
    long double inf = numeric_limits<double>::infinity();
    long double fact, i, res;
    fact = 1;
    for(i = 1; i < 10000.0; ++ i)
    {
        res = fact * i;
        cout << i << " " << setprecision(8) << res << endl;
        if (!(res < inf))
        {
            cout << endl << "inf" << endl;
            return 0;
        }
        fact = res;
    }
    return 0;
}
