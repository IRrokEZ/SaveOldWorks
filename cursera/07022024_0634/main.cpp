#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <set>

using namespace std;

int main() {
    map <set <string>, int> bus_stop;
    int Q, stop_count, counter = 1;
    string stop;
    set <string> stops;
    cin >> Q;
    for (int i = 0; i < Q; ++ i) {
        cin >> stop_count;
        stops.clear();
        for (int j = 0; j < stop_count; ++ j) {
            cin >> stop;
            stops.insert(stop);
        }
        if (bus_stop.find(stops) == bus_stop.end()) {
            bus_stop[stops] = counter;
            cout << "New bus " << counter << endl;
            ++ counter;
        } else {
            cout << "Already exists for " << bus_stop[stops] << endl;
        }
    }
    return 0;
}
