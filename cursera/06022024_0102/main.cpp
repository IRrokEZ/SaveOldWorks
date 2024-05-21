#include <iostream>
#include <string>
#include <vector>

using namespace std;

void MoveStrings (vector <string>& source, vector <string>& destination) {
    for (auto str : source) {
        destination.push_back(str);
    }
    source.clear();
}

int main()
{
    return 0;
}
