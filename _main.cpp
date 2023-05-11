#include <vector>
#include "matrix.h"
// #include "stats.h" //for stdz()
// #include <fstream>
#include <iostream>
// #include <iomanip>
// #include <map>
// #include <string>
// #include <utility>
// #include "rnt.h"
// #include "graph.h"
// #include "prob.h"
// #include "data.h"
#include "qr.h"
#include "qrc.h"
// #include <ctime>
// #include "lr.h"
// #include "svd.h"
// #include "svdc.h"
// #include <cmath>
// #include <vector>
// #include <string>
// #include "eig.h"
// #include "eigen.h"
#include "pca.h"
#include "data.h"
// #include "graph.h"
#include "TApplication.h"
using namespace std;

void StandaloneApplication(int argc, char **argv)
{
   Matrix A(3, 3, vector<double>{3, 14, 9, 6, 43, 3, 6, 22, 15});
   QR q(A);
   q.getQ().print();
   q.getR().print();
   (q.getQ() * q.getR()).print();
   cout << "\naspire" << endl;
}

int main(int argc, char **argv)
{

   TApplication app("ROOT Application", &argc, argv);
   StandaloneApplication(app.Argc(), app.Argv());
   app.Run();
   return 0;
}