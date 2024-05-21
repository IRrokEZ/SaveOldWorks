#include <iostream>
#include <string>

using namespace std;

bool IsPalindrom(string& s){
    if (s.size() < 2){
        return true;
    }
    for(size_t i = 0; i < s.size() / 2; ++ i) {
        if(s[i] != s[s.size() - i - 1]){
            return false;
        }
    }
    return true;
}

int main() {
    string s;
    cin >> s;
    cout << IsPalindrom(s);
    return 0;
}
