#include <chrono>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace std::chrono;

int main()
{
    duration<double> d{};
    cout << d.count() << endl;
    duration<int, milli> d1{1000};
    cout << fixed<<setprecision(10);
    cout << d1.count() << endl;
    return 0;
}
