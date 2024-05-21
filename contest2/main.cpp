#include <iostream>
#include <string>

using namespace std;

int getPalindrom(int number, size_t leading_zeros) {
    auto str = to_string(number);
    int result = 0;
    int pow = 1;
    for (int i = 0; i < str.length(); ++ i)
    {
        result += int(str.at(i) - '0') * pow;
        pow *= 10;
    }
    for (int i = 0; i < (leading_zeros - str.length()); ++ i)
    {
        result *= 10;
    }
    return result;
}

int main()
{
    int  mxHours, mxmns;
    cin >> mxHours >> mxmns;

    int mn, mx;
    if (mxmns > mxHours)
    {
        mn = mxHours - 1;
        mx = mxmns - 1;
    }
    else
    {
        mn = mxmns - 1;
        mx = mxHours - 1;
    }

    size_t leadingZeros = to_string(mx - 1).length();

    int c = 0;
    for (int i = 0; i <= mx; ++ i)
    {
        int pali = getPalindrom(i, leadingZeros);
        if (pali <= mn)
        {
            ++ c;
        }
    }
    cout << c << endl;
    return 0;
}
