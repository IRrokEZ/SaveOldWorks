#include <iostream>

using namespace std;

int Factorial(int a){
    if (a <= 1) {
        return 1;
    }
    int rez = 1;
    for(int i = 2; i <= a; ++ i) {
        rez *= i;
    }
    return rez;
}

int main() {
    int a;
    cin >> a;
    cout << Factorial(a);
    return 0;
}
