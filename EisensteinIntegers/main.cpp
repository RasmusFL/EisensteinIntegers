#include <iostream>
#include <windows.h>
#include "eisenstein_integers.h"
#include <time.h>

using namespace std;

int main()
{
    SetConsoleOutputCP(1253); // "ANSI" Greek
    EI alpha {2, 3};
    EI beta {-5, 8};
    cout << alpha + beta << endl;
    cout << alpha - beta << endl;
    cout << alpha * beta << endl;
    cout << beta / alpha << endl;
    cout << beta % alpha << endl;
    cout << EIgcd({36, -56}, {-144, 913}) << endl;
    vector<EI> v = EIexgcd({36, -56}, {-144, 913});
    for (EI a : v)
        cout << a << endl;
    cout << EImodlinearsolve({36, -56}, {-144, 913}, {587}) << endl;
    EImodlinearsolve({36, -56}, {-144, 913}, {8,2});
    printEIproduct(EIprimefactor({894582, -1870956}));
    cout << CubicResSymbol({36, -26}, {-144, 913}) << endl;
    return 0;
}
