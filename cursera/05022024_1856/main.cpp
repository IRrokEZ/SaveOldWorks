#include <iostream>
#include <string>

using namespace std;

int main(){
    string a, b, c;
    cin >> a >> b >> c;
    if((a <= b) && (a <= c)){
        cout << a;
        return 0;
    }
    if((b < a) && (b <= c)){
        cout << b;
        return 0;
    }
    if((c < a) && (c < b)){
        cout << c;
        return 0;
    }
    return 0;
}
