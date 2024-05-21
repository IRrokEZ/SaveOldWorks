#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

int main() {
    int a, b, c;
    cin >> a;
    if (a == 1) {
        cout << "1";
        return 0;
    }
    string rez = "";
    while (true) {
        b = a / 2;
        c = a % 2;
        if (c == 1){
            rez += "1";
        } else {
            rez += "0";
        }
        if (b == 1) {
            rez += "1";
            break;
        }
        a = b;
    }
    reverse(rez.begin(), rez.end());
    cout << rez;
    return 0;
}
