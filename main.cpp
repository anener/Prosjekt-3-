#include <iostream>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "time.h"


#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <omp.h>

#include "gauss_legendre.h"
#include "gauss_laguerre.h"
#include "montecarlo.h"

using namespace std;

void parallell(int);

int main()
{
    int n;
    double a, b;
    clock_t start, slutt;
    double tid;

    cout << "Antall integrasjons punkg (n): \n";
    cin >> n;

    a = -3.0;
    b = 3.0;

    //start = clock();
    //gaussLegendre(n, a, b);
    //slutt = clock();
    //tid = double(slutt - start)/CLOCKS_PER_SEC;
    //cout << "Tid = " << tid << "\n";

    //start = clock();
    //gaussLaguerre(n);
    //slutt = clock();
    //tid = double(slutt - start)/CLOCKS_PER_SEC;
    //cout << "Tid = " << tid << "\n";

    //start = clock();
    //uniform_moneCarlo(n, a, b);
    //slutt = clock();
    //tid = double(slutt - start)/CLOCKS_PER_SEC;
    //cout << "Tid = " << tid << "\n";

    start = clock();
    exp_monteCarlo(n);
    slutt = clock();
    tid = double(slutt - start)/CLOCKS_PER_SEC;
    cout << "Tid = " << tid << "\n";
    parallell(n);



}

//funksjon for aa parallellisere exp_monteCarlo()
void parallell(int n) {
    double wtime = omp_get_wtime(); //def tid

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

    int i;
    //starter parallelliseringen
    #pragma omp parallel for default(shared) private(i) reduction(+:mc)
    for(i=1; i<=n; i++) {
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

    wtime = omp_get_wtime() - wtime;
    cout << "Parallellisert tid = " << wtime << "\n";

    delete [] t; delete [] y; delete [] r; delete [] theta;
    delete [] phi;
}

// UNIT TESTING:
// intFunction

TEST_CASE("Tester intFunction") {
    //x1 = y1 = 2, z1 = 1
    //x2 = z2 = 2, y2 = 1
    double f = intFunction(2, 2, 1, 2, 1, 2);
    double x = exp(-24.0)/sqrt(2.0);

    REQUIRE(f == Approx(x));
}


TEST_CASE("Tester gaulag()") {
    int n = 20;
    double *x = new double [n];
    double *w = new double [n];
    double alpha = 3.0;
    gaulag(x, w, n, alpha);
    double f = 0.0;
    for(int i=1; i<= n; i++){
        f += w[i];
    }

    REQUIRE(f == Approx(6));
}


TEST_CASE("Tester gauleg()") {
    int n = 20;
    double pi = M_PI;

    double *phi = new double [n];
    double *w_phi = new double [n];
    gauleg(0.0, pi*2.0, phi, w_phi, n);

    double f = 0.0;
    for(int j=0; j<n; j++) {
        f += w_phi[j];
    }

    REQUIRE(f == Approx(2.0*pi));
}
