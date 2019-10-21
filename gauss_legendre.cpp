#include "gauss_legendre.h"

using namespace std;

//funksjon for aa regne ut Gauss-Legendre
void gaussLegendre(int n, double a, double b) {
    double *x = new double [n]; //lagre plass for abcissas verdier (mesh points)
    double *w = new double [n]; //lagre plass for weights verdier
    //a = start punkt for intervallet vi integrerer over
    //b = slutt punkt fonr intervallet vi integrerer over
    gauleg(a, b, x, w, n); //regner ut abcissas og weights av lengde n

    double gaussInt = 0.0;

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++) {
            for(int k=0; k<n; k++) {
                for(int l=0; l<n; l++) {
                    for(int p=0; p<n; p++) {
                        for(int q=0; q<n; q++) {
                            gaussInt += w[i]*w[j]*w[k]*w[l]*w[p]*w[q]*intFunction(x[i], x[j], x[k], x[l], x[p], x[q]); //for-loop for aa summere opp alle vekt verdiene og funksjonen med de forskjellige koordinatene
                        }
                    }
                }
            }
        }
    }
    cout << "integration = " << gaussInt << "\n";

    delete [] x;
    delete [] w;
}

//definerer funksjonene som skal integreres
double intFunction(double x1, double x2, double y1, double y2, double z1, double z2) {
    double alpha = 2;
    double r1, r2, exp1, exp2;
    r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
    r2 = sqrt(x2*x2 + y2*y2 + z2*z2);

    exp1 = -2*alpha*r1;
    exp2 = -2*alpha*r2;

    double deno = sqrt(pow((x1-x2), 2) + pow((y1-y2), 2) + pow((z1-z2), 2));
    double f;
    if(deno > 1e-10) {
        f = exp(exp1+exp2)/deno;
    }
    else {
        f = 0.0;
    }
    return f;
}
