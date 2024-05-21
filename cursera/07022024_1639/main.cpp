#include <iostream>
#include <vector>
#include <algorithm>
#include <cctype>

using namespace std;

void Print(const vector <pair <string, string>>& arr) {
    for (const auto& elt : arr) {
        cout << elt.second << " ";
    }
}

void Print2(const vector <pair <string, string>>& arr) {
    for (const auto& elt : arr) {
        cout << elt.first << " ";
    }
    cout << endl;
    for (const auto& elt : arr) {
        cout << elt.second << " ";
    }
    cout << endl;
}

int main() {
    int N;
    cin >> N;
    vector <pair <string, string>> ar(N);
    string str;
    for (int i = 0; i < N; ++ i) {
        cin >> str;
        ar[i].first = "";
        ar[i].second = str;
        for (const auto& el : ar[i].second) {
            ar[i].first += tolower(el);
        }
    }
    sort(ar.begin(), ar.end());
    Print(ar);
    return 0;
}
