#include <iostream>

using namespace std;

int main() {
    int a, b, c, d;
    cin >> a >> b;
    c = max(a, b);
    d = min(a, b);
    while (true) {
        int t = c % d;
        if (t == 0) {
            cout << d;
            break;
        } else {
            c = d;
            d = t;
        }
    }
    return 0;
}
