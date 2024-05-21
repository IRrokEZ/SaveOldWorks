#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main() {
    int N;
    cin >> N;
    if (N <= 0) {
        return 0;
    }
    vector <double> temp(N);
    string rez = "";
    double srar = 0;
    int counter = 0;
    for (int i = 0; i < N; ++ i) {
        cin >> temp[i];
        srar += temp[i];
    }
    srar /= N;
    for (int i = 0; i < N; ++ i) {
        if (temp[i] > srar) {
            ++ counter;
            rez += to_string(i) + " ";
        }
    }
    cout << counter << endl << rez;
    return 0;
}
