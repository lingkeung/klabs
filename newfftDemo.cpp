#include "newcmatrix.h"
#include "graph.h"
#include "TApplication.h"
#include "rnt.h"

void StandaloneApplication(int argc, char **argv)
{
    cout << "Spectrum Analysis Simulation" << endl;

    const double pi = M_PI;
    double fs = 500;               // sampling frequency [Hz]
    double Ns = pow(2, 14);        // no. of samples (radix-2)
    double T = Ns / fs;            // timespan [s]
    Matrix t = linspace(Ns, 0, T); // sampling time instants

    Matrix signal(Ns, 1); // test signal
    for (int i = 1; i <= Ns; i++)
    {
        signal(i) = sin(2 * pi * 50 * t(i)) + 0.5 * sin(2 * pi * 100 * t(i));
    }
    Random rand;
    Matrix noise = rand.matrix(Ns, 1, 'u', -5, 5); // random noise
    Matrix sPn = signal + noise;
    newcMatrix dft = (1 / Ns) * newfdft(sPn); // fdft() is a radix-2 fft function
    int N = Ns / 2 + 1;                 // display valid half of amplitude spectrum only
    Matrix spectrum(N, 1);              // prepare spectrum from dft
    for (int i = 1; i <= N; i++)
    {
        spectrum(i) = 2 * abs(dft(i));
    }

    double f1 = 1 / T;                             // fundamental frequency
    Matrix frequency = f1 * linspace(N, 0, N - 1); // prepare frequency range

    Graph gs(N, frequency, spectrum); // spectrum plot
    gs.set_opx("AL");                 // axes and line
    gs.set_title("Spectrum;Frequency [Hz];Amplitude [2*abs(Cn)]");
    gs.xy();

    Graph gt(Ns, t, sPn); // time amplitude plot
    gt.set_opx("AL");
    gt.set_title("Signal Plus Noise;Time [s];Amplitude");
    gt.xy();
}

int main(int argc, char **argv)
{
    TApplication app("ROOT Application", &argc, argv);
    StandaloneApplication(app.Argc(), app.Argv());
    app.Run();
    return 0;
}
