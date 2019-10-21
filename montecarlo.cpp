#include "montecarlo.h"
#include "gauss_legendre.h"
#include "gauss_laguerre.h"


//funksjon for aa regne ut Monte Carlo med uniform fordeling
void uniform_moneCarlo(int n, double a, double b) {
    double N = static_cast<double>(n);

    double *x = new double [6]; //(x, y, z) verdiene
    double *y = new double [n];

    double mc = 0.0;
    double sigma = 0.0;
    double f = 0.0;

    double jac = pow((b-a), 6);

    //finner tilfelldige verdier som er uniform fordelt
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> rr(0.0, 1.0);

    for(int i=1; i<=n; i++) {
        for(int j=0; j<6; j++) {
            x[j] = a + (b-a)*rr(gen);
        }
        f = intFunction(x[0], x[1], x[2], x[3], x[4], x[5]);
        mc += f;
        y[i] = f;
    }
    mc = mc/N;

    for(int i=1; i<=n; i++) {
        sigma += (y[i]-mc)*(y[i]-mc); // (x - mu)^2
    }
    sigma = sigma*jac/N;
    double std = sqrt(sigma)/sqrt(n);
    double integral = mc*jac;

    cout << "I = " << integral << "\n";
    cout << "std = " << std << "\n";
    delete [] x;
    delete [] y;
}


//finner Monte Carlo integral med eksponensiell fordeling
void exp_monteCarlo(int n) {
    double pi = M_PI;
    double N = static_cast<double>(n);

    double *t = new double [6];
    double *y = new double [n];

    double *r = new double [2];
    double *theta = new double [2];
    double *phi = new double [2];

    double f = 0.0;
    double mc = 0.0;
    double sigma = 0.0;
    double jac = (2.0*pi)*(2.0*pi)*pi*pi;

    //finner tilfelldige verdier som er uniform fordelt
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> rr(0.0, 1.0);

    for(int i=1; i<=n; i++) {
        for(int j=0; j<6; j++) {
            t[j] = rr(gen);
        }
        r[0] = -log(1.0-t[0]);
        r[1] = -log(1.0-t[1]);
        theta[0] = pi*t[2];
        theta[1] = pi*t[3];
        phi[0] = 2.0*pi*t[4];
        phi[1] = 2.0*pi*t[5];

        f = intFunction_exp(r[0], r[1], theta[0], theta[1], phi[0], phi[1]);
        mc += f;
        y[i] = f;
    }
    mc = mc/N;

    for(int i=1; i<=n; i++) {
        sigma += (y[i]-mc)*(y[i]-mc);
    }
    sigma = sigma*jac/N;
    double std = sqrt(sigma)/sqrt(n);

    double integral = mc*jac;
    cout << "I = " << integral << "\n";
    cout << "std = " << std << "\n";
}

//funksjon som skal integreres med eksponensiell fordeling
double intFunction_exp(double r1, double r2, double t1, double t2, double p1, double p2) {
    double alpha = 2.0;
    double cos_b = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
    double r12 = (r1*r1) + (r2*r2) - (2.0*r1*r2*cos_b);
    double dd = (r1*r1)*(r2*r2)*sin(t1)*sin(t2);
    double p = exp(-(r1+r2));
    double f;
    if(r12 > 1e-10) {
        f = (exp(-2.0*alpha*(r1+r2))*dd)/(p*sqrt(r12));
    }
    else {
        f = 0.0;
    }
    return f;
}
