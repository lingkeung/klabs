// #include <vector>
// #include "matrix.h"
// #include "stats.h" //for stdz()
// #include <fstream>
// #include <iostream>
// #include <iomanip>
// #include <map>
// #include <string>
// #include <utility>
// #include "rnt.h"
// #include "graph.h"
// #include "prob.h"
// #include "data.h"
// #include "qr.h"
// #include "qrc.h"
// #include <ctime>
// #include "lr.h"
// #include "svd.h"
// #include "svdc.h"
// #include <cmath>
// #include <vector>
// #include <string>
// #include "eig.h"
// #include "eigen.h"
// #include "pca.h"
#include "data.h"
// #include "graph.h"
// #include "cn.h"
#include "TApplication.h"

using namespace std;

void StandaloneApplication(int argc, char **argv)
{
   DataR dr("water.dat", "waterR.dat");
   for (auto &&i : dr.a)
   {
      cout << i << "  " << endl;
   }
   for (auto &&i : dr.r)
   {
      cout << i << "  " << endl;
   }
}

int main(int argc, char **argv)
{
   TApplication app("ROOT Application", &argc, argv);
   StandaloneApplication(app.Argc(), app.Argv());
   app.Run();
   return 0;
}