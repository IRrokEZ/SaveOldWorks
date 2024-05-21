#include <iostream>
#include <cstdlib>
#include <time.h>

using namespace std;

int main()
{
    srand(time(NULL));
    for(int i = 0; i < 20; i ++)
    {
        int a = rand() % 10 - 5;
        cout << a << endl;
    }
    return 0;
}
