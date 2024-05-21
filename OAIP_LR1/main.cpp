#include <iostream>
#include <clocale>
#include <math.h>
#include <limits>
#include <iomanip>

using namespace std;

long double Factorial (int base)
{
    if (base > 1754)
    {
        cout << endl << "  ���������� ���������� ����� " << base << " ����������. ������������ ����������." << endl;
        return -1;
    }

    if ((base == 0) || (base == 1))
    {
        return 1;
    }
    long double result = 1;
    for (long double i = 1; i < (base + 1); ++ i)
    {
        result = result * i;
    }
    return result;
}

long double Func (long double x, int k, long double fact)
{
    long double koefB, convertK, result;
    if ((k % 2) == 0)
    {
        koefB = -1;
    }
    else
    {
        koefB = fact;
    }
    convertK = k;
    result = sin(convertK * x) / convertK * koefB;
    return result;
}

int main()
{
    setlocale(LC_ALL,"Russian");
    cout << endl << "  ����������, ������� �������� ������ ���������." << endl;
    long double x = 0, eps = -1, fact, valueFunk, previousvalueFunc, delta;
    long double inf = numeric_limits<double>::infinity();
    long double Summ = 0;
    int a = 12000000;
    int64_t termcounter = 0;
    while ((!(x > 0)) && (!(x < 0)))
    {
        cout << "  ������� �" << endl << "  ";
        cin >> x;
        if ((!(x > 0)) && (!(x < 0)))
        {
            cout << endl << "  ������� ������������ ��������. � �� ������ ���� ����� 0." << endl;
        }
    }
    while (!(eps > 0))
    {
        cout << "  ������� Epsilon" << endl << "  ";
        cin >> eps;
        if (!(eps > 0))
        {
            cout << endl << "  ������� ������������ ��������. Epsilon ������ ���� ������ ������ 0." << endl;
        }
    }
    while (abs(a) > 1000000)
    {
        cout << "  ������� �" << endl << "  ";
        cin >> a;
        if (abs(a) > 1000000)
        {
            cout << endl << "  ������� ������������ ��������. ����� � �� ������ �� ������ ������������ 1000000 (���� �������)." << endl;
        }
    }
    fact = Factorial(a);
    if (fact < 0)
    {
        return 0;
    }
    previousvalueFunc = eps;
    valueFunk = eps;
    delta = 2 * eps;
    for (int64_t i = 1; (abs(delta) > eps); ++ i)
    {
        long double tempFuncvalue = Func(x, i, fact);
        if ((tempFuncvalue < inf) && ((tempFuncvalue > 0) || (tempFuncvalue < 0)))
        {
            previousvalueFunc = valueFunk;
            valueFunk = tempFuncvalue;
            Summ += valueFunk;
            termcounter ++;
            delta = previousvalueFunc - valueFunk;
            if (abs(delta) < eps)
            {
                break;
            }
        }
    }
    cout << endl << "  �����:" << endl << "  ��� �������� �������� eps = " << fixed << showpoint << setprecision(8) << eps << "," << endl;;
    cout << "  �������� � = " << fixed << showpoint << setprecision(8) << x << "," << endl << "  �������� a = " << a << endl;
    cout << "  C���� ���� ����� " << fixed << showpoint << setprecision(8) << Summ << endl;
    cout << "  ���-�� ���������: " << termcounter;

    return 0;
}
