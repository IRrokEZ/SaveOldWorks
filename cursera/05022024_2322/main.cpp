#include <iostream>
#include <string>
#include <vector>

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

vector <string> PalindromFilter(vector <string> &words, int minLength){
    size_t sz;
    vector <string> rez = {};
    if (minLength >= 0) {
        sz = minLength;
    } else {
        sz = 0;
    }
    for (auto str : words) {
        if (str.size() >= sz) {
            if (IsPalindrom(str)) {
                rez.push_back(str);
            }
        }
    }
    return rez;
}


int main()
{

    return 0;
}
