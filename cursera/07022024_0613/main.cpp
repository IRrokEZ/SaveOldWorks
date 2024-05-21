#include <iostream>
#include <set>
#include <string>
#include <map>

using namespace std;

set<string> BuildMapValuesSet(const map<int, string>& m) {
    set <string> result;
    for (const auto& pair : m) {
        result.insert(pair.second);
    }
    return result;
}
