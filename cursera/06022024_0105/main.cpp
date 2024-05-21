#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

void Reverse(vector<int>& v) {
    size_t vsize = v.size();
    for(size_t i = 0; i < vsize / 2; ++ i) {
        swap(v[i], v[vsize - i - 1]);
    }
}

int main()
{
    return 0;
}
