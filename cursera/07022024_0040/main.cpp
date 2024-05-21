#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <set>

using namespace std;

int main() {
    map <string, vector <string>> bus_stop; // bus + stop_array
    map <string, vector <string>> stop_bus; // stop + all buses
    int Q, stop_count;
    string command, bus, stop_name, last_bus;
    cin >> Q;
    for (int i = 0; i < Q; ++ i) {
        cin >> command;
        if (command == "NEW_BUS") {
            cin >> bus >> stop_count;
            last_bus = bus;
            for (int j = 0; j < stop_count; ++ j) {
                cin >> stop_name;
                bus_stop[bus].push_back(stop_name);
                stop_bus[stop_name].push_back(bus);
            }
        } else if (command == "BUSES_FOR_STOP") {
            cin >> stop_name;
            if (stop_bus.find(stop_name) == stop_bus.end()) {
                cout << "No stop" << endl;
            } else {
                // cout << "Stop " << stop_name << ": ";
                for (const auto& st : stop_bus[stop_name]) {
                    cout << st << " ";
                }
                cout << endl;
            }
        } else if (command == "STOPS_FOR_BUS") {
            cin >> bus;
            if (bus_stop.find(bus) == bus_stop.end()) {
                cout << "No bus" << endl;
            } else {
                for (const auto& bs : bus_stop[bus]) {
                    cout << "Stop " + bs << ": ";
                    if (stop_bus[bs].size() == 1) {
                        cout << "no interchange" << endl;
                    } else {
                        for (const auto& st: stop_bus[bs]) {
                            if (st != bus) {
                                cout << st << " ";
                            }
                        }
                        cout << endl;
                    }
                }
            }
        } else if (command == "ALL_BUSES") {
            if (bus_stop.empty()) {
                cout << "No buses" << endl;
            } else {
                for (const auto& pair : bus_stop) {
                    cout << "Bus " << pair.first << ": ";
                    for (const auto& stp : pair.second) {
                        cout << stp << " ";
                    }
                    cout << endl;
                }
            }
        }
    }
    return 0;
}
