#include <iostream>
#include <vector>

using namespace std;

vector<int> Reversed(const vector<int>& v) {
    vector <int> rez = {};
    if (v.size() <= 1) {
        rez = v;
        return rez;
    }
    for (size_t i = v.size(); i > 0; --i) {
        rez.push_back(v[i - 1]);
    }
    return rez;
}

int main()
{

    return 0;
}
