#include <iostream>
#include <vector>
#include <algorithm>
#include <cctype>
#include <set>

using namespace std;

class SortedStrings {
public:
    void AddString(const string& s) {
        my_strings.insert(s);
    }
    vector<string> GetSortedStrings() {
        vector <string> result (my_strings.size());
        int cnt = 0;
        for (const auto& str : my_strings) {
            result[cnt] = str;
            ++ cnt;
        }
        return result;
  }
private:
    multiset <string> my_strings;
};

void Print(const vector <string>& arr) {
    for (const auto& elt : arr) {
        cout << elt << " ";
    }
    cout << endl;
}

int main() {
    SortedStrings strings;
    strings.AddString("first");
    strings.AddString("third");
    strings.AddString("second");
    Print(strings.GetSortedStrings());
    strings.AddString("second");
    Print(strings.GetSortedStrings());
    return 0;
}
