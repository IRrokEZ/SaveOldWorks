#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

void Print(const vector <int>& arr) {
    for (const auto& elt : arr) {
        cout << elt << " ";
    }
}

int main() {
    int N;
    cin >> N;
    vector <int> ar(N);
    for (auto &el : ar) {
        cin >> el;
    }
    sort(ar.begin(), ar.end(), [](int& a, int& b){return abs(a) < abs(b);});
    Print(ar);
    return 0;
}
