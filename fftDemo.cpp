#include "cmatrix.h"
#include "graph.h"
#include "TApplication.h"

void StandaloneApplication(int argc, char **argv)
{
    cout << "Spectrum Analysis Simulation" << endl;

    const double pi = M_PI;
    double fs = 500;              // sampling frequency
    double N = pow(2, 14);        // no. of samples (radix-2)
    double T = N / fs;            // timespan
    Matrix t = linspace(N, 0, T); // sampling instants

    Matrix signal(N, 1); // test signal
    for (int i = 1; i <= N; i++)
    {
        signal(i) = sin(2 * pi * 50 * t(i)) + 0.5 * sin(2 * pi * 100 * t(i));
    }
    cMatrix dft = (1 / N) * fdft(signal); // fdft() is a radix-2 fft function

    N = N / 2 + 1; // display valid half of spectrum only

    Matrix spectrum(N, 1); // prepare spectrum
    for (int i = 1; i <= N; i++)
    {
        spectrum(i) = 2 * abs(dft(i));
    }

    double f1 = 1 / T;                             // fundamental frequency
    Matrix frequency = f1 * linspace(N, 0, N - 1); // prepare frequency range

    Graph g(N, frequency, spectrum); // plot spectrum
    g.xy();
}

int main(int argc, char **argv)
{
    TApplication app("ROOT Application", &argc, argv);
    StandaloneApplication(app.Argc(), app.Argv());
    app.Run();
    return 0;
}
