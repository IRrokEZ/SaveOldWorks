#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main() {
    vector <int> que;
    int Q, i = 0, arg;
    string command;
    cin >> Q;
    while (i < Q) {
        cin >> command;
        if (command == "WORRY") {
            cin >> arg;
            que[arg] = 1;
        } else if (command == "QUIET") {
            cin >> arg;
            que[arg] = 0;
        } else if (command == "COME") {
            cin >> arg;
            if (arg > 0) {
                for (int p = 0; p < arg; ++ p) {
                    que.push_back(0);
                }
            } else {
                arg *= -1;
                for (int p = 0; p < arg; ++ p) {
                    que.pop_back();
                }
            }
        } else if (command == "WORRY_COUNT") {
            int counter = 0;
            for (auto p : que) {
                if (p == 1) {
                    ++ counter;
                }
            }
            cout << counter << endl;
        }
        ++ i;
    }
    return 0;
}
