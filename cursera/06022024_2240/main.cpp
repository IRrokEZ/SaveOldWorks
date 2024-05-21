#include <iostream>
#include <map>
#include <vector>
#include <string>

using namespace std;

int main() {
    map <string, string> my_map;
    int Q;
    string command, p1, p2, save;
    vector <string> result;
    cin >> Q;
    for(int i = 0; i < Q; ++ i) {
        cin >> command;
        if (command == "CHANGE_CAPITAL") {
            cin >> p1 >> p2;
            if (my_map.find(p1) == my_map.end()) {
                my_map[p1] = p2;
                result.push_back("Introduce new country " + p1 + " with capital " + p2);
            } else {
                if (my_map[p1] == p2) {
                    result.push_back("Country " + p1 + " hasn't changed its capital");
                } else {
                    result.push_back("Country  " + p1 + "  has changed its capital from " + my_map[p1]
                                     + " to " + p2);
                    my_map[p1] = p2;
                }
            }
        } else if (command == "RENAME") {
            cin >> p1 >> p2;
            if ((p1 == p2)
               || (my_map.find(p2) != my_map.end())
               || (my_map.find(p1) == my_map.end())) {
                result.push_back("Incorrect rename, skip");
            } else {
                save = my_map[p1];
                my_map.erase(p1);
                my_map[p2] = save;
                result.push_back("Country " + p1 + " with capital " + save + " has been renamed to " + p2);
            }
        } else if (command == "ABOUT") {
            cin >> p1;
            if (my_map.find(p1) == my_map.end()) {
                result.push_back("Country " + p1 + " doesn't exist");
            } else {
                result.push_back("Country " + p1 + " has capital " + my_map[p1]);
            }
        } else if (command == "DUMP") {
            if (my_map.empty()) {
                result.push_back("There are no countries in the world");
            } else {
                for (const auto& pair : my_map) {
                    result.push_back(pair.first + "/" + pair.second);
                }
            }
        }
    }
    for (const auto& s : result) {
        cout << s << endl;
    }
    return 0;
}
