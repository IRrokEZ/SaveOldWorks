#include <iostream>
#include <vector>
#include <string>
#include <utility>

using namespace std;

const vector <size_t> days_in_month = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; // шаблон кол-ва дней в месяце

int main() {
    int Q, i = 0, current_month = 0;
    vector <pair <int, vector <string>>> timesheet (31, make_pair(0, vector <string> ())); // вектор пар чисо, список дел
    string command, name;
    size_t arg, days_before, dif, current_days;

    for (size_t j = 0; j < days_in_month[current_month]; ++ j) {
        timesheet[j].first = j + 1;
    }

    cin >> Q;

    while (i < Q) {
        cin >> command;
        if (command == "ADD") {
            cin >> arg >> name;
            timesheet[arg - 1].second.push_back(name);
        } else if (command == "DUMP") {
            cin >> arg;
            cout << timesheet[arg - 1].second.size() << " ";
            for (const auto& wrk : timesheet[arg - 1].second) {
                cout << wrk << " ";
            }
            cout << endl;
        } else if (command == "NEXT") {
            days_before = days_in_month[current_month];
            ++ current_month;
            if (current_month > 11) {
                current_month = 0;
            }
            current_days = days_in_month[current_month];
            if (current_days < days_before) {
                dif = days_before - current_days;
                for (size_t p = 0; p < dif; ++ p) {
                    timesheet[current_days - 1].second.insert(end(timesheet[current_days - 1].second),
                                                              begin(timesheet[current_days + p].second),
                                                              end(timesheet[current_days + p].second));
                    timesheet[current_days + p].second.clear();
                }
            }
        }
        ++ i;
    }
    return 0;
}
