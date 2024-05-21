#include <iostream>
#include <string>

using namespace std;

int main() {
    string inp;
    cin >> inp;
    int fcount = 0;
    for (int i = 0; i < inp.size(); ++ i) {
        if (inp[i] == 'f') {
            ++ fcount;
        }
        if (fcount == 2) {
            cout << i;
            return 0;
        }
    }
    if (fcount == 1) {
        cout << -1;
        return 0;
    }
    cout << -2;
    return 0;
}
