#include <iostream>
#include <cmath>

using namespace std;

int main(){
    double a, b, c;
    cin >> a >> b >> c;
    if (a != 0) {
        double D = b * b - 4 * a * c;
        if (D > 0) {
            cout << (-b + sqrt(D)) / 2.0 / a << " " << (-b - sqrt(D)) / 2.0 / a;
        } else if (D < 0) {
            return 0;
        } else {
            cout << - (b / 2.0 / a);
        }
    } else {
        if (b != 0) {
            cout << -(c / b);
        }
    }
    return 0;
}
