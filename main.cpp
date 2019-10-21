#include <iostream>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "time.h"

#include <omp.h>

#include "gauss_legendre.h"
#include "gauss_laguerre.h"
#include "montecarlo.h"

using namespace std;

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

    //start = clock();
    //exp_monteCarlo(n);
    //slutt = clock();
    //tid = double(slutt - start)/CLOCKS_PER_SEC;
    //cout << "Tid = " << tid << "\n";

    int thread_num = omp_get_max_threads();
        cout << "Number of processors: " << omp_get_num_procs() << "\n";
        cout << "Number og threads: " << thread_num << "\n";

        double wtime = omp_get_wtime(); //def tid
        double *x = new double [n];
        double *w = new double [n];
        gauleg(a, b, x, w, n);
        double gaussInt = 0.0;
    #pragma ompi parallel for default(shared) private(i, j, k, l, p, q) reduction(+:gaussInt)
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++) {
                for(int k=0; k<n; k++) {
                    for(int l=0; l<n; l++) {
                        for(int p=0; p<n; p++) {
                            for(int q=0; q<n; q++) {
                                gaussInt += w[i]*w[j]*w[k]*w[l]*w[p]*w[q]*intFunction(x[i], x[j], x[k], x[l], x[p], x[q]);
                            }
                        }
                    }
                }
            }
        }
        wtime = omp_get_wtime() - wtime;
        cout << wtime << "\n";
        cout << gaussInt << "\n";



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
