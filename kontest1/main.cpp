#include <iostream>
#include <vector>

using namespace std;

int main()
{
    int n, q;
    cin >> n;
    if ((n % 2) == 1)
    {
        for(int i = 0; i < n; i ++)
        {
            cin >> q;
        }
        cout << -1;
        return 0;
    }
    vector <int> a(n / 2);
    for(int i = 0; i < n / 2; i ++)
    {
        cin >> a[i];
    }
    for(int i = n / 2 - 1; i >= 0; -- i)
    {
        cin >> q;
        a[i] += q;
    }
    q = 0;
    for(int i = 1; i < n/2; i ++)
    {
        if(a[i - 1] != a[i])
        {
            cout << -1;
            return 0;
        }
    }
    cout << a[0];
    return 0;
}
