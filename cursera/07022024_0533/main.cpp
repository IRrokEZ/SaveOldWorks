#include <iostream>
#include <set>
#include <string>

using namespace std;

int main() {
    set <string> ms;
    string t;
    int n;
    cin >> n;
    for (int i = 0; i < n; ++ i) {
        cin >> t;
        ms.insert(t);
    }
    cout << ms.size();
    return 0;
}
