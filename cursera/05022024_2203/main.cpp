#include <iostream>

using namespace std;

int main() {
    double n, a, b, x, y;
    cin >> n >> a >> b >> x >> y;
    if (n > b) {
        cout << n - n * y / 100.0;
        return 0;
    }
    if (n > a) {
        cout << n - n * x / 100.0;
        return 0;
    }
    cout << n;
    return 0;
}
