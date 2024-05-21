#include <iostream>
#include <set>
#include <string>
#include <map>

using namespace std;

int main() {
    int Q;
    string command, word1, word2;
    map <string, set <string>> dict;
    cin >> Q;
    for (int i = 0; i < Q; ++ i) {
        cin >> command;
        if (command == "ADD") {
            cin >> word1 >> word2;
            dict[word1].insert(word2);
            dict[word2].insert(word1);
        } else if (command == "COUNT") {
            cin >> word1;
            cout << dict[word1].size() << endl;
        } else if (command == "CHECK") {
            cin >> word1 >> word2;
            bool flag = true;
            if (dict[word1].find(word2) == dict[word1].end()) {
                flag = false;
            } else if (dict[word2].find(word1) == dict[word2].end()) {
                flag = false;
            }
            if (flag) {
                cout << "YES" << endl;
            } else {
                cout << "NO" << endl;
            }
        }
    }
    return 0;
}
