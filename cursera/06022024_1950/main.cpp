#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

int main() {
    map <char, int> letters;
    int N, let_counter;
    bool flag;
    string first, second;
    vector <string> answer;
    cin >> N;
    for (int i = 0; i < N; ++ i) {
        flag = true;
        letters.clear();
        let_counter = 0;
        cin >> first >> second;
        for (size_t j = 0; j < first.size(); ++ j) {
            ++ letters[first[j]];
            ++ let_counter;
        }
        for (size_t j = 0; j < second.size(); ++ j) {
            -- letters[second[j]];
            if (letters[second[j]] < 0) {
                flag = false;
            }
            -- let_counter;
        }
        if ((flag) && (let_counter == 0)) {
            answer.push_back("YES");
        } else {
            answer.push_back("NO");
        }
    }
    for (const auto& el : answer) {
        cout << el << endl;
    }
    return 0;
}
