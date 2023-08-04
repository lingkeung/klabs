#include <iostream>
#include <iomanip>

using namespace std;

class X
{
public:
    X(){};
    X(const X &)
    {
        cout << "Copy Constructor function" << endl;
    }
    X &operator=(const X &)
    {
        cout << "Assignment copy function" << endl;
        return *this;
    }
};

X fun()
{
    X t;
    return t;
}

int main()
{
    X x;
    // X y = x;
    x = fun();
}
